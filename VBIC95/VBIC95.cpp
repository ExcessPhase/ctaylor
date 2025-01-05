#ifdef __JACOBIAN__
#include "../cjacobian.h"
#else
#include "../ctaylor.h"
#endif
#include "../LUFAC/lufac.h"
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <cassert>
#define SELF_HEATING
namespace vbic95
{
using namespace boost::mp11;
enum class enumParameters
{
#define __create__(a, b) a
#define __create2__(a, b) a##_TNOM
#define __COMMA__ ,
#include "parameters.h"
};
enum class enumMembers
{
#define __create__(a) a
#define __COMMA__ ,
#include "members.h"
};
/// for indexing of const char* key inside std::map
struct compare
{	bool operator()(const char *const _p0, const char *const _p1) const
	{	return std::strcmp(_p0, _p1) < 0;
	}
};
/// translates parameter name into enumeration
/// fails with an exception
static enumParameters getParameterIdByName(const char *const _p)
{	static const std::map<const char*, enumParameters, compare> s = {
#define __create__(a, b) {#a, enumParameters::a}
#define __create2__(a, b) {#a, enumParameters::a##_TNOM}
#define __COMMA__ ,
#include "parameters.h"
	};
	return s.at(_p);
}
#ifdef __JACOBIAN__
using namespace jacobian;
#else
using namespace taylor;
#endif
/// part of the original VBIC standard
template<typename P_T, typename EA_T, typename VTV_T, typename RT_T>
auto psibi (const P_T& P, const EA_T&EA, const VTV_T&Vtv, const RT_T&rT) 
{
	const auto psiio = 2.0 * Vtv * log( exp ( 0.5 * P / Vtv ) - exp ( - 0.5 * P / Vtv ) );
	const auto psiin = psiio * rT - 3.0 * Vtv * log ( rT ) - EA * ( rT - 1.0 );
	return psiin + 2.0 * Vtv * log ( 0.5 * ( 1.0 + sqrt ( 1.0 + 4.0 * exp ( - psiin / Vtv ) ) ) );
}
/// part of the original VBIC standard
template<typename V_T, typename P_T, typename M_T, typename FC_T, typename A_T>
auto qj(const V_T&V, const P_T&P, const M_T&M, const FC_T&FC, const A_T&A )
{
	return if_(
		A <= 0.0,
		[&](void)
		{
			//
			//		SPICE regional depletion capacitance model
			//
			const auto dvh = V - FC * P;
			return if_(
				dvh > 0.0,
				[&](void)
				{	const auto qlo = P * ( 1.0 - pow( 1.0 - FC, 1.0 - M ) ) / ( 1.0 - M );
					const auto qhi = dvh * ( 1.0 - FC + 0.5 * M * dvh / P ) / ( pow( 1.0 - FC, 1.0 + M ) );
					return qlo + qhi;
				},
				[&](void)
				{	const auto qlo = P * ( 1.0 - pow( 1.0 - V / P, 1.0 - M ) ) / ( 1.0 - M );
					const auto qhi = 0.0;
					return qlo + qhi;
				}
			);
		},
		[&](void)
		{
	//
	//		Single piece depletion capacitance model
	//
	//		Based on c=1/(1-V/P)^M, with sqrt limiting to make it
	//		C-inf continuous (and computationally efficient), with
	//		added terms to make it monotonically increasing (which
	//		is physically incorrect, but avoids numerical problems
	//		and kinks near where the depletion and diffusion
	//		capacitances are of similar magnitude), and with appropriate
	//		offsets added to that qj(V=0)=0.
	//
			const auto dv0 = - P * FC;
			const auto mv0 = sqrt ( dv0 * dv0 + A );
			const auto vl0 = 0.5 * ( dv0 - mv0 ) + P * FC;
			const auto q0 = - P * pow( 1.0 - vl0 / P, 1.0 - M ) / ( 1.0 - M );
			const auto dv = V - P * FC;
			const auto mv = sqrt ( dv * dv + A );
			const auto vl = 0.5 * ( dv - mv ) + P * FC;
			const auto qlo = - P * pow( 1.0 - vl / P, 1.0 - M ) / ( 1.0 - M );
			return qlo + pow( 1.0 - FC, - M ) * ( V - vl + vl0 ) - q0;
		}
	);
}
/// part of the original VBIC standard
template<typename V_T, typename P_T, typename M_T, typename AV1_T, typename AV2_T>
auto avalm(const V_T&V, const P_T&P, const M_T&M, const AV1_T&AV1, const AV2_T&AV2)
{
//
//		Kloosterman/de Graaff weak avalanche model
//
	const auto vl = 0.5 * ( sqrt ( sqr( P - V ) + 0.01 ) + ( P - V ) );
	return AV1 * vl * exp ( - AV2 * pow(vl, M - 1.0 ) );
}
typedef std::map<enumMembers, double> enumMembers2double;
static constexpr double KB = 1.380662E-23;	// Boltzmann's constant, J/K
static constexpr double QQ = 1.602189E-19;	// magnitude of electronic charge, C
static constexpr double TABS = 2.731500E+02;	// 0 degrees C in degrees K
static constexpr double TAMB = 27.0;
/// reads parameters from an external file
/// called from the constructor of VBIC
static enumMembers2double readParams(const char *const _p)
{	std::ifstream sFile(_p);
	sFile.exceptions(std::ios_base::badbit | std::ios_base::failbit);
	sFile.exceptions(std::ios_base::goodbit);
	double d;
	std::string s;
	enumMembers2double sRet;
	std::map<enumParameters, double> sParameters;
#define __create__(a, b) sParameters.insert(std::make_pair(enumParameters::a, b));
#define __create2__(a, b) sParameters.insert(std::make_pair(enumParameters::a##_TNOM, b));
#define __COMMA__
#include "parameters.h"
	while (true)
	{	if (sFile >> d >> s)
			sParameters[getParameterIdByName(s.c_str())] = d;
		else
			break;
	}

#ifndef SELF_HEATING
#define __create__(a, b) const auto a = sParameters.at(enumParameters::a);
#define __create2__(a, b) const auto a##_TNOM = sParameters.at(enumParameters::a##_TNOM);
#define __COMMA__
#include "parameters.h"
#define __create__(a, b) const auto a = b;
#define __create2__(a, b) const auto a = b;
#define __COMMA__
#include "temperatureSetup.h"
#define __create__(a) sRet[enumMembers::a] = a;
#define __COMMA__
#include "members.h"
#else
#define __create__(a) sRet[enumMembers::a] = sParameters.at(enumParameters::a);
#define __COMMA__
#include "members.h"
#endif
	return sRet;
}
/// enumeration for circuit nodes
/// negative nodes are not part of the matrix and have a nonzero voltage associated with them
enum class enumCircuitNodes:std::ptrdiff_t
{
#define __create__(a) a
#define __create2__(a, b) a = b
#define __COMMA__ ,
#include "circuitNodes.h"
};
/// nodes which are part of the transistor
enum class enumNodes
{
#define __create__(a) a
#define __COMMA__ ,
#include "nodes.h"
};
/// names of transistor nodes
static const char*const s_aNodeNames[] =
{
#define __create__(a) #a
#define __COMMA__ ,
#include "nodes.h"
};
/// names of circuit nodes
static const char*const s_aCircuitNodeNames[] =
{
#define __create__(a) #a,
#define __create2__(a, b)
#define __COMMA__
#include "circuitNodes.h"
};
/// input voltages branches
enum class enumBranches:std::size_t
{
#define __COMMA__ ,
#define __create__(a, b, c) a
#include "inputs.h"
};
/// input voltage branches as pairs of transistor nodes
static constexpr const std::pair<enumNodes, enumNodes> s_aInput2NodePair[] =
{
#define __COMMA__ ,
#define __create__(a, b, c) {enumNodes::b, enumNodes::c}
#include "inputs.h"
};
/// current sources inside the transistor
enum class enumCurrentOutputs:std::size_t
{
#define __create__(a, b, c) a
#define __COMMA__ ,
#include "currentSources.h"
};
/// are we using cjacobian.h or ctaylor.h
/// definition of the input voltage type
#ifdef __JACOBIAN__
template<enumBranches E>
using createIndep = cjacobian<
	mp_list<
		mp_size_t<std::size_t(E)>
	>
>;
#else
	/// !
	/// here one needs to change to 2 in order to use Halley's method
	/// which is incompatible with the delta-x scaling
constexpr std::size_t MAX = 2;
template<enumBranches E>
using createIndep = ctaylor<makeIndependent<std::size_t(E)>, MAX>;
#endif
/// the transistor instance
struct vbic
{
#define __create__(a) const double a;
#define __COMMA__
#include "members.h"
		/// the constructor
	vbic(const enumMembers2double &_r)
		:
#define __create__(a) a(_r.at(enumMembers::a))
#define __COMMA__ ,
#include "members.h"
	{
	}
		/// access of the input voltage of a transistor node
	template<enumNodes _e>
	inline static double getInputVoltage(
			/// fixed voltages for certain nodes not part of the circuit matrix
		const std::array<double, 2>&_rN,
			/// the node voltages part of the circuit matrix
		const std::array<double, std::ptrdiff_t(enumCircuitNodes::NumberOfNodes)>& _rP,
			/// mapping of internal nodes to external ones
			/// negative enums are connected to fixed voltages
		const std::array<std::ptrdiff_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> & _pT
	)
	{	if (_e == enumNodes::NumberOfNodes)
			return 0.0;
		else
		{	const auto i = _pT.at(std::size_t(_e));
			if (i < 0)
				return _rN.at(-1 - i);
			else
			if (i == std::ptrdiff_t(enumCircuitNodes::NumberOfNodes))
				return 0.0;
			else
				return _rP.at(i);
		}
	}
		/// the same for node pairs
	template<enumNodes _e0, enumNodes _e1>
	inline static double getInputVoltage(
		const std::array<double, 2>&_rN,
		const std::array<double, std::ptrdiff_t(enumCircuitNodes::NumberOfNodes)>& _rP,
		const std::array<std::ptrdiff_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> & _pT
	)
	{	return getInputVoltage<_e0>(_rN, _rP, _pT) - getInputVoltage<_e1>(_rN, _rP, _pT);
	}
		/// the entry point for calculation
	auto calculate(
			/// node voltages part of the matrix
		const std::array<double, std::ptrdiff_t(enumCircuitNodes::NumberOfNodes)>& _rNodeVoltages,
			/// fixed node voltages
		const std::array<double, std::size_t(2)>& _rNodeVoltages2,
			/// the 3-dimensional hessian
		lufac::index2Index2Index2Double&_r2,
			/// the 2-dimensional jacobian
		lufac::index2Index2Double&_rO,
			/// the 1-dim output value
		lufac::index2Double &_rV,
			/// translation of the internal transistor nodes to circuit Nodes
		const std::array<std::ptrdiff_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> & _pT,
			/// for returning the DC terminal currents
		std::array<double, 4>&_rTC
	) const
	{
		/// creating the input voltages with correct type
#define __create__(Vbe, b, e) const auto Vbe = createIndep<enumBranches::Vbe>(\
		getInputVoltage<enumNodes::b, enumNodes::e>(_rNodeVoltages2, _rNodeVoltages, _pT),\
		false\
	);
#define __COMMA__
#include "inputs.h"
		/// self-heating code
#ifdef SELF_HEATING
#define __create__(a, b) const auto a = b;
#define __create2__(a, b) const auto a = b;
#define __COMMA__
#include "temperatureSetup.h"
#endif

//	This section defines branch currents, charges. and fluxes as functions
//	of the branch voltages, currents, and model parameters.

//	First some bias-independent calculations

	const auto IVEF = VEF <= 0.0 ? 0.0 : 1.0 / VEF;
	const auto IVER = VER <= 0.0 ? 0.0 : 1.0 / VER;
	const auto IIKF = IKF <= 0.0 ? 0.0 : 1.0 / IKF;
	const auto IIKR = IKR <= 0.0 ? 0.0 : 1.0 / IKR;
	const auto IIKP = IKP <= 0.0 ? 0.0 : 1.0 / IKP;
	const auto IVO = VO <= 0.0 ? 0.0 : 1.0 / VO;
	const auto IHRCF = HRCF <= 0.0 ? 0.0 : 1.0 / HRCF;
	const auto IVTF = VTF <= 0.0 ? 0.0 : 1.0 / VTF;
	const auto IITF = ITF <= 0.0 ? 0.0 : 1.0/ITF;
	const auto slTF = ITF <= 0.0 ? 1.0 : 0.0;
#ifdef EXCESS_PHASE
	const auto LEP = TD <= 0.0 ? 0.0 : TD/3.0;
	const auto CEP = TD <= 0.0 ? 0.0 : TD;
#endif

//	Calculate normalized depletion charges

	const auto qdbe  = qj ( Vbei, PE, ME, FC, AJE );	// b-e depletion charge
	const auto qdbex = qj ( Vbex, PE, ME, FC, AJE );	// b-e depletion charge (side)
	const auto qdbc  = qj ( Vbci, PC, MC, FC, AJC );	// b-c depletion charge
	const auto qdbep = qj ( Vbep, PC, MC, FC, AJC );	// parasitic b-e deplt'n charge
					// Note that b-c and parasitic b-e
					// junctions are taken to have same
					// built-in potential and grading coeff.
	const auto qdbcp = qj ( Vbcp, PS, MS, FC, AJS );	// parasitic b-c deplt'n charge

//	Transport currents and qb

	const auto Itfi  = IS * ( exp ( Vbei / ( NF * Vtv ) ) - 1.0 );	// fwd ideal
	const auto Itri  = IS * ( exp ( Vbci / ( NR * Vtv ) ) - 1.0 );	// rev ideal
	const auto q1z   = 1.0 + qdbe * IVER + qdbc * IVEF;			// proper q1
	const auto q1    = 0.5 * ( sqrt ( sqr( q1z - 0.0001 ) + 0.0001 * 0.0001 )
			+ q1z - 0.0001 ) + 0.0001;		// 1e-4 limit
	const auto q2    = Itfi * IIKF + Itri * IIKR;
	const auto qb    = 0.5 * ( q1 + sqrt ( q1 * q1 + 4.0 * q2 ) );	// no limiting
	const auto Itzr  = Itri / qb;					// static rev
	const auto Itzf  = Itfi / qb;					// static fwd
#ifdef EXCESS_PHASE
	const auto Flxf  = LEP * Ilxf;					// flux in Lxf
	const auto Qcxf  = CEP * Vcxf;					// charge in Cxf
	const auto Itxf  = Vrxf;						// ex-ph fwd
#endif

//	Transport currents and qb of the parasitic, to give Iccp

	const auto Itfp  = ISP * ( WSP * exp ( Vbep / ( NFP * Vtv ) )
			+ ( 1.0 - WSP ) * exp ( Vbci / ( NFP * Vtv ) ) - 1.0 );
	assert(isfinite(Itfp) && !isnan(Itfp));
	const auto Itrp  = ISP * ( exp ( Vbcp / ( NFP * Vtv ) ) - 1.0 );
	assert(isfinite(Itrp) && !isnan(Itrp));
	const auto q2p   = Itfp * IIKP;				// only fwd part
	assert(isfinite(q2p) && !isnan(q2p));
	const auto qbp   = 0.5 * ( 1.0 + sqrt ( 1.0 + 4.0 * q2p ) );// no Early effect
	assert(isfinite(qbp) && !isnan(qbp));
	const auto Iccp  = ( Itfp - Itrp ) / qbp;			// parasitic transport I
	assert(isfinite(Iccp) && !isnan(Iccp));

//	Currents in resistors, with nodes collapsed for zero resistance
//	Note that RBI, RCI and RBP are modulated

	const auto Ircx = if_(
		RCX <= 0.0,
		[](void)
		{	std::abort();
			return 0.0;
		},
		[&](void)
		{	return Vrcx / RCX;
		}
	);
	const auto Irci_Kbci_Kbcx = if_(
		RCI <= 0.0,
		[](void)
		{	std::abort();
			return std::make_tuple(0.0, 0.0, 0.0);
		},
		[&](void)
		{	const auto Vbcx = Vbci - Vrci;
			assert(isfinite(Vbcx) && !isnan(Vbcx));
			const auto Kbci = sqrt ( 1.0 + GAMM * exp ( Vbci / Vtv ) );
			assert(isfinite(Kbci) && !isnan(Kbci));
			const auto Kbcx = sqrt ( 1.0 + GAMM * exp ( Vbcx / Vtv ) );
			assert(isfinite(Kbcx) && !isnan(Kbcx));
			const auto rKp1 = ( Kbci + 1.0 ) / ( Kbcx + 1.0 );
			assert(isfinite(rKp1) && !isnan(rKp1));
			const auto Iohm = ( Vrci + Vtv * ( Kbci - Kbcx - log ( rKp1 ) ) ) / RCI;
			assert(isfinite(Iohm) && !isnan(Iohm));
			const auto derf = IVO * RCI * Iohm / ( 1.0 + 0.5 * IVO * IHRCF *
				sqrt ( Vrci * Vrci +  0.01 ) );
			assert(isfinite(derf) && !isnan(derf));
			const auto Irci = Iohm / sqrt ( 1.0 + derf * derf  );
			assert(isfinite(Irci) && !isnan(Irci));
			return std::make_tuple(Irci, Kbci, Kbcx);
		}
	);
	const auto &Irci = std::get<0>(Irci_Kbci_Kbcx);
	const auto &Kbci = std::get<1>(Irci_Kbci_Kbcx);
	const auto &Kbcx = std::get<2>(Irci_Kbci_Kbcx);
	const auto Irbx = if_(
		RBX  <= 0.0,
		[](void)
		{	std::abort();
			return 0.0;
		},
		[&](void)
		{	return Vrbx / RBX;
		}
	);
	const auto Irbi = if_(
		RBI  <= 0.0,
		[](void)
		{	std::abort();
			return 0.0;
		},
		[&](void)
		{	return Vrbi * qb / RBI;	// simple qb modulation model for now,
					// other models to be defined.
					// Irbi=Irbi(Vrbi,Vbei,Vbci) because
					// qb  =qb  (     Vbei,Vbci)
		}
	);
	const auto Ire = if_(
		RE   <= 0.0,
		[](void)
		{	std::abort();
			return 0.0;
		},
		[&](void)
		{	return Vre  / RE;
		}
	);
	const auto Irs = if_(
		RS   <= 0.0,
		[](void)
		{	std::abort();
			return 0.0;
		},
		[&](void)
		{	return Vrs  / RS;
		}
	);
	const auto Irbp = if_(
		RBP <= 0.0,
		[](void)
		{	std::abort();
			return 0.0;
		},
		[&](void)
		{	return Vrbp * qbp / RBP;	// Irbp=Irbp(Vrbp,Vbep) because
						// qbp = qbp(     Vbep)
		}
	);

//	b-e and b-c components of base current of intrinsic device

	const auto Ibe   = WBE *
		(   IBEI  * ( exp ( Vbei / ( NEI  * Vtv ) ) - 1.0 )
		  + IBEN  * ( exp ( Vbei / ( NEN  * Vtv ) ) - 1.0 ) );
	const auto Ibex  = ( 1.0 - WBE ) *
		(   IBEI  * ( exp ( Vbex / ( NEI  * Vtv ) ) - 1.0 )
		  + IBEN  * ( exp ( Vbex / ( NEN  * Vtv ) ) - 1.0 ) );
	const auto Ibc   = (   IBCI  * ( exp ( Vbci / ( NCI  * Vtv ) ) - 1.0 )
		  + IBCN  * ( exp ( Vbci / ( NCN  * Vtv ) ) - 1.0 ) );

//	Parasitic b-e current, with calculation bypass for efficiency
//	Emission coefficients are same as for intrinsic b-c, as they
//	are the same junction
	const auto Ibep = if_(
		IBEIP <= 0.0 && IBENP <= 0.0,
		[](void)
		{	return 0.0;
		},
		[&](void)
		{	return IBEIP * ( exp ( Vbep / ( NCI * Vtv ) ) - 1.0 )
				+ IBENP * ( exp ( Vbep / ( NCN * Vtv ) ) - 1.0 );
		}
	);

//	Parasitic b-c current, with calculation bypass for efficiency.
//	This element should normally never be of importance, but is
//	included to detect incorrect biasing, a useful task.
	const auto Ibcp = if_(
		IBCIP <= 0.0 && IBCNP <= 0.0,
		[](void)
		{	return 0.0;
		},
		[&](void)
		{	return IBCIP * ( exp ( Vbcp / ( NCIP * Vtv ) ) - 1.0 )
				+ IBCNP * ( exp ( Vbcp / ( NCNP * Vtv ) ) - 1.0 );
		}
	);

//	b-c weak avalanche current
	const auto Igc = if_(
		AVC1 <= 0.0,
		[](void)
		{	return 0.0;
		},
		[&](void)
		{	return ( Itzf - Itzr - Ibc ) * avalm ( Vbci, PC, MC, AVC1, AVC2 );
		}
	);

//	Charge elements, Qbe diffusion charge is not split at present

//	This is the only place where C-inf continuity is broken, as the
//	SPICE forward transit time bias dependence, which includes "if"
//	conditions, is used
	const auto sgItf = Itfi > 0.0 ?  1.0 : 0.0;
	const auto rItf  = Itfi * sgItf * IITF;
	const auto tff   = TF * ( 1.0 + QTF * q1 ) *
		( 1.0 + XTF * exp ( Vbci * IVTF / 1.44 )
			* ( slTF + sqr( rItf / ( rItf + 1.0 ) ) ) * sgItf );
	const auto Qbe   = CJE  * qdbe  * WBE           + tff * Itfi / qb;
	const auto Qbex  = CJE  * qdbex * ( 1.0 - WBE );
	const auto Qbc   = CJC  * qdbc                  + TR  * Itri      + QCO * Kbci;
	const auto Qbcx  =                                                  QCO * Kbcx;
	const auto Qbep  = CJEP * qdbep                 + TR  * Itfp;
	const auto Qbcp  = CJCP * qdbcp;			// no diffusion charge
	const auto Qbeo  = CBEO * Vbe;			// extrinsic b-e overlap charge
	const auto Qbco  = CBCO * Vbc;			// extrinsic b-c overlap charge

#ifdef SELF_HEATING

//	Thermal power generation must be done by summing I*V over all
//	non-energy storage elements of the electrical model

	const auto Ith_Irth_Qcth = if_(
		RTH  <= 0.0,
		[](void)
		{	return std::make_tuple(0.0, 0.0, 0.0);
		},
		[&](void)
		{	const auto Ith   =   max(Ibe  * Vbei, 0.0) + max(Ibc  * Vbci, 0.0)
			+ max(( Itzf - Itzr ) * ( Vbei - Vbci ), 0.0)
			+ max(Ibep * Vbep, 0.0) + max(Ibcp * Vbcp, 0.0) + max(Iccp * ( Vbep - Vbcp ), 0.0)
			+ max(Ircx * Vrcx, 0.0) + max(Irci * Vrci, 0.0) + max(Irbx * Vrbx, 0.0)
			+ max(Irbi * Vrbi, 0.0) + max(Ire  * Vre, 0.0)  + max(Irbp * Vrbp, 0.0)
			+ max(Irs  * Vrs, 0.0)  + max(Ibex * Vbex, 0.0) + max(- Igc  * Vbci, 0.0);

//		Simple linear thermal resistance and capacitance, could be
//		made nonlinear if necessary

			const auto Irth  = delT / RTH;

			const auto Qcth  = CTH  * delT;
			return std::make_tuple(Ith, Irth, Qcth);
		}
	);
	const auto &Ith = std::get<0>(Ith_Irth_Qcth);
	const auto &Irth = std::get<1>(Ith_Irth_Qcth);
	const auto &Qcth = std::get<2>(Ith_Irth_Qcth);
#endif
		/// creating code to write into RHS and jacobian and hessian
#define __create__(a, b, c)\
{	{	const auto iBT = _pT[std::size_t(enumNodes::b)];\
		if (iBT >= 0 && iBT != std::ptrdiff_t(enumCircuitNodes::NumberOfNodes))\
			writeOutput(_rV[iBT], _rO[iBT], _r2[iBT], -a, _pT);\
	}\
	{	const auto iCT = _pT[std::size_t(enumNodes::c)];\
		if (iCT >= 0 && iCT != std::ptrdiff_t(enumCircuitNodes::NumberOfNodes))\
			writeOutput(_rV[iCT], _rO[iCT], _r2[iCT], a, _pT);\
	}\
}
#define __COMMA__
#include "currentSources.h"
//	End of equations.
			/// code for calculating DC terminal currents
		_rTC[0]=value(Itzf-Itzr-Ibc+Igc+Irbp);
		_rTC[1]=value(Ibe+Ibex+Ibc-Igc+Ibep+Iccp);
		_rTC[2]=value(-Itzf+Itzr-Ibe-Ibex);
		_rTC[3]=-_rTC[0]-_rTC[1]-_rTC[2];
	}
		/// instantiated in writeOutput with the template argument being identical to the argument to cjacobian/ctaylor
	template<typename T>
	struct copyOutput
	{	lufac::index2Double&m_rO;
		lufac::index2Index2Double&m_r1;
		double&m_rValue;
#ifdef __JACOBIAN__
		const cjacobian<T>&m_rV;
#else
		const ctaylor<T, MAX>&m_rV;
#endif
		const std::array<std::ptrdiff_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> &m_pT;
			/// constructor
		copyOutput(
			double&_rValue,
			lufac::index2Double&_rO,
			lufac::index2Index2Double&_r1,
#ifdef __JACOBIAN__
			const cjacobian<T>&_rV,
#else
			const ctaylor<T, MAX>&_rV,
#endif
			const std::array<std::ptrdiff_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> &_pT
		)
			:m_rO(_rO),
			m_r1(_r1),
			m_rV(_rV),
			m_rValue(_rValue),
			m_pT(_pT)
		{
#ifdef __JACOBIAN__
				/// writing the RHS which is not part of the template argument in case of jacobian.h
			m_rValue += m_rV.m_s.back();
#endif
		}
			/// translating internal nodes to external nodes
		constexpr enumCircuitNodes translate(const enumNodes _e) const
		{	const auto i = m_pT[std::size_t(_e)];
			if (i < 0)
				return enumCircuitNodes::NumberOfNodes;
			else
			{	assert(i <= std::ptrdiff_t(enumCircuitNodes::NumberOfNodes));
				assert(i != std::ptrdiff_t(18446744073709551614ul));
				return enumCircuitNodes(i);
			}
		}
			/// dito for a pair of nodes
		constexpr std::pair<enumCircuitNodes, enumCircuitNodes> translate(const std::pair<enumNodes, enumNodes>&_r) const
		{	return std::make_pair(translate(_r.first), translate(_r.second));
		}
			/// callback for second order derivatives
		template<std::size_t POS>
		void operator()(const enumCircuitNodes _e0, const enumCircuitNodes _e1, const mp_size_t<POS>&, const bool _bMinus, const bool _bDouble = false) const
		{	assert(std::ptrdiff_t(_e0) >= std::ptrdiff_t());
			assert(std::ptrdiff_t(_e1) >= std::ptrdiff_t());
			if (_e0 != enumCircuitNodes::NumberOfNodes && _e1 != enumCircuitNodes::NumberOfNodes)
				if (const double d = m_rV.m_s.at(POS))
				{	auto &r = m_r1[std::size_t(_e0)][std::size_t(_e1)];
					if (_bDouble)
						if (_bMinus)
							r -= 2.0*d;
						else
							r += 2.0*d;
					else
						if (_bMinus)
							r -= d;
						else
							r += d;
				}
		}
#ifndef __JACOBIAN__
			/// second order same derivative
		template<typename ENUM>
		void operator()(const mp_list<pair<ENUM, mp_size_t<2> > >&) const
		{	const auto sNodePair = translate(s_aInput2NodePair[ENUM::value]);
			typedef mp_list<pair<ENUM, mp_size_t<1> > > LIST;
			constexpr std::size_t POS = mp_find<T, LIST>::value;
			(*this)(sNodePair.first, sNodePair.first, mp_size_t<POS>(), false, true);
			(*this)(sNodePair.second, sNodePair.second, mp_size_t<POS>(), false, true);
			(*this)(sNodePair.first, sNodePair.second, mp_size_t<POS>(), true, true);
			(*this)(sNodePair.second, sNodePair.first, mp_size_t<POS>(), true, true);
		}
			/// second order cross derivatives
		template<typename ENUM0, typename ENUM1>
		void operator()(const mp_list<pair<ENUM0, mp_size_t<1> >, pair<ENUM1, mp_size_t<1> > >&) const
		{	const auto sNodePair0 = translate(s_aInput2NodePair[ENUM0::value]);
			const auto sNodePair1 = translate(s_aInput2NodePair[ENUM1::value]);
			typedef mp_list<pair<ENUM0, mp_size_t<1> >, pair<ENUM1, mp_size_t<1> > > LIST;
			constexpr std::size_t POS = mp_find<T, LIST>::value;

			(*this)(sNodePair0.first, sNodePair1.first, mp_size_t<POS>(), false);
			(*this)(sNodePair0.second, sNodePair1.second, mp_size_t<POS>(), false);
			(*this)(sNodePair0.first, sNodePair1.second, mp_size_t<POS>(), true);
			(*this)(sNodePair0.second, sNodePair1.first, mp_size_t<POS>(), true);
			
			(*this)(sNodePair1.first, sNodePair0.first, mp_size_t<POS>(), false);
			(*this)(sNodePair1.second, sNodePair0.second, mp_size_t<POS>(), false);
			(*this)(sNodePair1.first, sNodePair0.second, mp_size_t<POS>(), true);
			(*this)(sNodePair1.second, sNodePair0.first, mp_size_t<POS>(), true);
		}
			/// first order
		template<typename ENUM>
		void operator()(const mp_list<pair<ENUM, mp_size_t<1> > >&) const
		{	const auto sNodePair = translate(s_aInput2NodePair[ENUM::value]);
			typedef mp_list<pair<ENUM, mp_size_t<1> > > LIST;
			constexpr std::size_t POS = mp_find<T, LIST>::value;
			if (const double d = m_rV.m_s.at(POS))
			{	if (sNodePair.first != enumCircuitNodes::NumberOfNodes)
					m_rO[std::size_t(sNodePair.first)] -= d;
				if (sNodePair.second != enumCircuitNodes::NumberOfNodes)
					m_rO[std::size_t(sNodePair.second)] += d;
			}
		}
			/// 0th order
		void operator()(const mp_list<>&) const
		{	m_rValue += m_rV.m_s.at(mp_find<T, mp_list<> >::value);
		}
			/// in case of someone sets MAX to something larger than 2
		template<typename ...PAIRS>
		void operator()(const mp_list<PAIRS...>&) const
		{	static_assert(implementation::order<mp_list<PAIRS...> >::value > 2, "order must be larger than 2!");
		}
#else
			/// for jacobian.h
		template<std::size_t ENUM>
		void operator()(const mp_size_t<ENUM>&) const
		{	const auto sNodePair = translate(s_aInput2NodePair[ENUM]);
			constexpr std::size_t POS = mp_find<T, mp_size_t<ENUM> >::value;
			if (const double d = m_rV.m_s.at(POS))
			{	if (sNodePair.first != enumCircuitNodes::NumberOfNodes)
					m_rO[std::size_t(sNodePair.first)] -= d;
				if (sNodePair.second != enumCircuitNodes::NumberOfNodes)
					m_rO[std::size_t(sNodePair.second)] += d;
			}
		}
#endif
	};
		/// writing output
	template<typename T>
	static void writeOutput(
		double&_rValue,
		lufac::index2Double&_rO,
		lufac::index2Index2Double&_r1,
#ifdef __JACOBIAN__
		const cjacobian<T>&_rV,
#else
		const ctaylor<T, MAX>&_rV,
#endif
		const std::array<std::ptrdiff_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> & _pT
	)
	{	assert(isfinite(_rV) && !isnan(_rV));
		assert(std::isfinite(_rValue) && !std::isnan(_rValue));
		mp_for_each<T>(copyOutput<T>(_rValue, _rO, _r1, _rV, _pT));
	}
	static const char*const s_aNames[];
};
const char*const vbic::s_aNames[] = {
#define __create__(a, b, c) #a
#define __COMMA__ ,
#include "currentSources.h"
};
thread_local std::size_t s_iIndent = 0;
std::ostream&printIdent(std::ostream&_rS)
{	for (std::size_t i = 0; i < s_iIndent; ++i)
		_rS << "        ";
	return _rS;
}
template<typename ...REST>
std::ostream &print(std::ostream &_rS, const std::tuple<REST...>&, const std::integral_constant<std::size_t, 0>&)
{	return _rS;
}
template<std::size_t MPOS, typename ...REST>
std::ostream &print(std::ostream &_rS, const std::tuple<REST...>&_r, const std::integral_constant<std::size_t, MPOS>&)
{	_rS << vbic::s_aNames[sizeof...(REST) - MPOS] << "=" << std::get<sizeof...(REST) - MPOS>(_r) << ",";
	return print(_rS, _r, std::integral_constant<std::size_t, MPOS-1>());
}
template<typename ...REST>
std::ostream &operator<<(std::ostream &_rS, const std::tuple<REST...>&_r)
{	return print(_rS, _r, std::integral_constant<std::size_t, sizeof...(REST)>());
}
std::ostream &operator<<(std::ostream &_rS, const enumNodes _e)
{	return _rS << s_aNodeNames[std::size_t(_e)];
}
template<typename T0, typename T1>
std::ostream &operator<<(std::ostream &_rS, const std::pair<T0, T1>&_r)
{	return printIdent(_rS) << _r.first << "," << _r.second << "\n";
}
template<typename K, typename V>
std::ostream &operator<<(std::ostream &_rS, const std::map<K, V>&_r)
{	++s_iIndent;
	for (const auto &r : _r)
		_rS << r << "\n";
	--s_iIndent;
	return _rS;
}
/// translation of internal transistor nodes to external circuit nodes
constexpr enumCircuitNodes translateNodes(const enumNodes _e)
{	switch (_e)
	{	default:
			throw std::logic_error("Invalid node ID!");
		case enumNodes::NumberOfNodes:
			return enumCircuitNodes::NumberOfNodes;
		case enumNodes::c:
			return enumCircuitNodes::c;
		case enumNodes::b:
			return enumCircuitNodes::b;
		case enumNodes::e:
			return enumCircuitNodes::NumberOfNodes;
		case enumNodes::s:
			return enumCircuitNodes::NumberOfNodes;
#ifdef SELF_HEATING
		case enumNodes::dt:
			return enumCircuitNodes::dt;
		case enumNodes::tl:
			return enumCircuitNodes::NumberOfNodes;
#endif
		case enumNodes::cx:
			return enumCircuitNodes::cx;
		case enumNodes::ci:
			return enumCircuitNodes::ci;
		case enumNodes::bx:
			return enumCircuitNodes::bx;
		case enumNodes::bi:
			return enumCircuitNodes::bi;
		case enumNodes::ei:
			return enumCircuitNodes::ei;
		case enumNodes::si:
			return enumCircuitNodes::si;
		case enumNodes::bp:
			return enumCircuitNodes::bp;
#ifdef EXCESS_PHASE
		case enumNodes::xf1:
			return enumCircuitNodes::xf1;
		case enumNodes::xf2:
			return enumCircuitNodes::xf2;
#endif
	}
}
lufac::index2Double operator+(const lufac::index2Double&_r0, const lufac::index2Double&_r1)
{	lufac::index2Double s;
	lufac::index2Double::const_iterator pR0, pR1;
	for (pR0 = _r0.begin(), pR1 = _r1.begin(); pR0 != _r0.end() && pR1 != _r1.end(); )
		if (pR0->first < pR1->first)
		{	s[pR0->first] = pR0->second;
			++pR0;
		}
		else
		if (pR0->first > pR1->first)
		{	s[pR1->first] = pR1->second;
			++pR1;
		}
		else
		{	s[pR1->first] = pR0->second + pR1->second;
			++pR0;
			++pR1;
		}
	for (; pR0 != _r0.end(); ++pR0)
		s[pR0->first] = pR0->second;
	for (; pR1 != _r1.end(); ++pR1)
		s[pR1->first] = pR1->second;
	return s;
}
lufac::index2Index2Double operator+(const lufac::index2Index2Double&_r0, const lufac::index2Index2Double&_r1)
{	lufac::index2Index2Double sRet;
	lufac::index2Index2Double::const_iterator pR0, pR1;
	for (pR0 = _r0.begin(), pR1 = _r1.begin(); pR0 != _r0.end() && pR1 != _r1.end(); )
		if (pR0->first < pR1->first)
		{	sRet[pR0->first] = pR0->second;
			++pR0;
		}
		else
		if (pR0->first > pR1->first)
		{	sRet[pR1->first] = pR1->second;
			++pR1;
		}
		else
		{	sRet[pR0->first] = pR0->second + pR1->second;
			++pR0;
			++pR1;
		}
	for (; pR0 != _r0.end(); ++pR0)
		sRet[pR0->first] = pR0->second;
	for (; pR1 != _r1.end(); ++pR1)
		sRet[pR1->first] = pR1->second;
	return sRet;
}
/// transposed(_rL)*_rR
lufac::index2Double operator*(const double _d, const lufac::index2Double&_r)
{	lufac::index2Double s;
	for (const auto &r : _r)
		s[r.first] = r.second*_d;
	return s;
}
lufac::index2Index2Double operator*(const lufac::index2Index2Index2Double&_rL, const lufac::index2Double&_rR)
{	lufac::index2Index2Double s;
	for (const auto &r0 : _rL)
	{	auto &rS0 = s[r0.first];
		for (const auto &r1 : r0.second)
		{	auto &rS1 = rS0[r1.first];
			for (const auto &r2 : r1.second)
			{	const auto pFind = _rR.find(r2.first);
				if (pFind != _rR.end())
					rS1 += r2.second*pFind->second;
			}
		}
	}
	return s;
}
}
int main(int argc, char**argv)
try
{	using namespace vbic95;
	using namespace lufac;
	if (argc < 2)
	{	std::cerr << argv[0] << ": usage: " << argv[0] << " PARS" << std::endl;
		return 1;
	}
		/// single transistor instance
	const vbic sI(readParams(argv[1]));
		/// translation of internal nodes to external ones
	const std::array<std::ptrdiff_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> sTrans{
#define __create__(a) std::ptrdiff_t(translateNodes(enumNodes::a))
#define __COMMA__ ,
#include "nodes.h"
	};
		/// the voltage sweep
	for (double vb = 0.7; vb <= 1.00001; vb += 0.05)
		for (double vc = 0.0; vc <= 5.00001; vc += 0.05)
		{		/// the external voltages (RHS)
			std::array<double, std::ptrdiff_t(enumCircuitNodes::NumberOfNodes)> sV({});
				/// guess for internal voltages to be identical to external ones
			sV[std::ptrdiff_t(enumCircuitNodes::bx)] = vb;
			sV[std::ptrdiff_t(enumCircuitNodes::bi)] = vb;
			sV[std::ptrdiff_t(enumCircuitNodes::bp)] = vc;
			sV[std::ptrdiff_t(enumCircuitNodes::ei)] = 0.0;
			sV[std::ptrdiff_t(enumCircuitNodes::si)] = 0.0;
			sV[std::ptrdiff_t(enumCircuitNodes::ci)] = vc;
			sV[std::ptrdiff_t(enumCircuitNodes::cx)] = vc;
	//ve=0;vs=0
			//for(vb=0.7;vb<=1.00001;vb+=0.05)
			//for(vc=0.0;vc<=5.00001;vc+=0.05)
		//print vc,vb,ve,vs
				/// terminal currents
			std::array<double, 4> sTC;
			for (std::size_t i = 0; i < 200; ++i)
			{	index2Index2Index2Double sH;
				index2Index2Double sJ;
				index2Double sValues;
				static_assert(std::ptrdiff_t(enumCircuitNodes::b) == -1, "circuit node b must be negative!");
				static_assert(std::ptrdiff_t(enumCircuitNodes::c) == -2, "circuit node c must be negative!");
					/// calling the single transistor instance
					/// second argument are the fixed external voltages (with negative indicies)
				sI.calculate(sV, {vb, vc}, sH, sJ, sValues, sTrans, sTC);
				for (const auto &r0 : sH)
					for (const auto &r1 : r0.second)
						for (const auto &r2 : r1.second)
							if (!std::isfinite(r2.second) || std::isnan(r2.second))
								throw std::logic_error("Calculate produces NAN!");
				for (const auto &r0 : sJ)
					for (const auto &r1 : r0.second)
						if (!std::isfinite(r1.second) || std::isnan(r1.second))
							throw std::logic_error("Calculate produces NAN!");
				for (const auto &r : sValues)
					if (!std::isfinite(r.second) || std::isnan(r.second))
						throw std::logic_error("Calculate produces NAN!");
/*
x_{k+1} = x_k - [J(x_k) - (1/2) H(x_k) F(x_k)]^{-1} F(x_k)
*/
				if (sValues.size() != sJ.size())
					throw std::logic_error("matrix size is not identical!");
				const auto sJ1 = sJ + sH*(0.5*sValues);
				const auto sFactored = factor(sJ1, sValues);
				const auto sDelta = solve(sFactored, sValues);
				for (const auto &r : sDelta)
					if (!std::isfinite(r.second) || std::isnan(r.second))
						throw std::logic_error("solve produces NAN!");
					/// calculating the norm of the RHS
					/// not divided by number to be identical to original solver.f
				const auto dNormO = std::sqrt(
					std::accumulate(
						sValues.begin(),
						sValues.end(),
						0.0,
						[](const double _dSum, const index2Double::value_type&_r)
						{	return _dSum + _r.second*_r.second;
						}
					)
				);
#if 0
				const double vscale = 1.0;
#else
				const auto vscale = [&](void)
				{	double vscale = 1.0;
					const auto dvmax = 0.2;
					for (const auto &r : sDelta)
#ifdef SELF_HEATING
						if (r.first != std::size_t(enumCircuitNodes::dt))
#endif
						if (std::abs(r.second)/dvmax*vscale > 1.0)
							vscale = std::abs(dvmax/r.second);
					return vscale;
				}();
#endif
					/// calculating the norm of the resulting delta
					/// not divided by number to be identical to original solver.f
				const auto dNorm = std::sqrt(
					std::accumulate(
						sDelta.begin(),
						sDelta.end(),
						0.0,
						[vscale](const double _dSum, const index2Double::value_type&_r)
						{
#ifdef SELF_HEATING
							if (_r.first == std::size_t(enumCircuitNodes::dt))
								return _dSum;
#endif
							const double d = _r.second*vscale;
							return _dSum + d*d;
						}
					)
				);
					/// applying the delta
					/// making certain temperature can not become negative
					/// the transistor does not cool the environment
				for (const auto &r : sDelta)
				{	double &rD = sV[r.first] += vscale*r.second;
#ifdef SELF_HEATING
					if (r.first == std::size_t(enumCircuitNodes::dt))
						if (rD < 0.0)
							rD = 0.0;
#endif
				}
					/// are we converged
				if (dNorm < 1e-6 && dNormO < 1e-6)
					break;
			}
				/// display the result in the same way as the original solver.f
#ifdef SELF_HEATING
			std::cout << vc << "\t" << vb << "\t" << sTC[0] << "\t" << sTC[1] << "\t" << -sV[std::size_t(enumCircuitNodes::dt)];
#else
			std::cout << vc << "\t" << vb << "\t" << sTC[0] << "\t" << sTC[1] << "\t" << sTC[3];
#endif
			std::cout << std::endl;
#ifdef SELF_HEATING
	//write(6,*) 'Inputs  are: v(c) v(b) v(e) v(s)'
	//write(6,*) 'Outputs are: i(c) i(b) i(e) dt'
//Vce Vbe Ic Ib Is
#else
	//write(6,*) 'Inputs  are: v(c) v(b) v(e) v(s)'
	//write(6,*) 'Outputs are: i(c) i(b) i(e) i(s)'
//Vce Vbe Ic Ib Is
#endif
		}
} catch (const std::exception&_r)
{	std::cerr << argv[1] << ": Exception caught: " << _r.what() << std::endl;
	return 1;
}
