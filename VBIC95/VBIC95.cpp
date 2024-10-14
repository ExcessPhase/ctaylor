#include "../ctaylor.h"
#include "../LUFAC/lufac.h"
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <numeric>
#include <algorithm>
//#define SELF_HEATING
//#define __DIODE__
namespace vbic95
{
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
struct compare
{	bool operator()(const char *const _p0, const char *const _p1) const
	{	return std::strcmp(_p0, _p1) < 0;
	}
};
static enumParameters getParameterIdByName(const char *const _p)
{	static const std::map<const char*, enumParameters, compare> s = {
#define __create__(a, b) {#a, enumParameters::a}
#define __create2__(a, b) {#a, enumParameters::a##_TNOM}
#define __COMMA__ ,
#include "parameters.h"
	};
	return s.at(_p);
}
using namespace taylor;
template<typename P_T, typename EA_T, typename VTV_T, typename RT_T>
auto psibi (const P_T& P, const EA_T&EA, const VTV_T&Vtv, const RT_T&rT) 
{
	const auto psiio = 2.0 * Vtv * log( exp ( 0.5 * P / Vtv ) - exp ( - 0.5 * P / Vtv ) );
	const auto psiin = psiio * rT - 3.0 * Vtv * log ( rT ) - EA * ( rT - 1.0 );
	return psiin + 2.0 * Vtv * log ( 0.5 * ( 1.0 + sqrt ( 1.0 + 4.0 * exp ( - psiin / Vtv ) ) ) );
}
template<typename V_T, typename P_T, typename M_T, typename FC_T, typename A_T>
auto qj(const V_T&V, const P_T&P, const M_T&M, const FC_T&FC, const A_T&A )
{
	return taylor::if_(
		A <= 0.0,
		[&](void)
		{
			//
			//		SPICE regional depletion capacitance model
			//
			const auto dvh = V - FC * P;
			return taylor::if_(
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
static const double LOGSQRT2 = std::log(std::sqrt(2.0));
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

enum class enumCircuitNodes
{
#define __create__(a) a
#define __COMMA__ ,
#include "circuitNodes.h"
};
enum class enumNodes
{
#define __create__(a) a
#define __COMMA__ ,
#include "nodes.h"
};
static const char*const s_aNodeNames[] =
{
#define __create__(a) #a
#define __COMMA__ ,
#include "nodes.h"
};
enum class enumBranches:std::size_t
{
#define __COMMA__ ,
#define __create__(a, b, c) a
#include "inputs.h"
};
static constexpr const std::pair<enumNodes, enumNodes> s_aInput2NodePair[] =
{
#define __COMMA__ ,
#define __create__(a, b, c) {enumNodes::b, enumNodes::c}
#include "inputs.h"
};
enum class enumCurrentOutputs:std::size_t
{
#define __create__(a, b, c) a
#define __COMMA__ ,
#include "currentSources.h"
};
template<enumBranches E>
using createIndep = ctaylor<makeIndependent<std::size_t(E)>, 2>;
struct vbic
{
#define __create__(a) const double a;
#define __COMMA__
#include "members.h"
	vbic(const enumMembers2double &_r)
		:
#define __create__(a) a(_r.at(enumMembers::a))
#define __COMMA__ ,
#include "members.h"
	{
	}
#if 0
	template<typename T>
	static auto EXP(const ctaylor<T, 1>&vnew, const double vold, bool *const _pbLimit) -> decltype(exp(vnew))
	{	if (!_pbLimit)
			return exp(vnew);
		else
#if 1
		if (_pbLimit && vnew > vold + 1.0e-1)
		{	*_pbLimit = true;
			const auto YOLD = std::exp(vold);
			const auto YPOLD = YOLD;
			const auto YINTENDED = YOLD + YPOLD*(value(vnew) - vold);
			const auto XINTENDED = log(YINTENDED);
			return YINTENDED*(vnew - value(vnew) + 1.0) - YINTENDED*(XINTENDED - value(vnew));
		}
		else
			return exp(vnew);
#else
		{	const auto vcrit = log(1.0/(std::sqrt(2.0)*IS));
			if (vnew > vcrit && abs(vnew - vold) > 2.0)
			{	*_pbLimit = true;
				if (vold > 0.0)
				{	const auto arg = 1.0 + vnew - vold;
					if (arg > 0)
						return exp(vold + log(arg));
					else
						return exp(vcrit);
				}
				else
					return vnew;
			}
			else
				return exp(vnew);
		}
#if 0
	if (vnew > vcrit && std::abs(vnew - vold) > vt + vt)
	{	*icheck = 1;
		if (vold > 0)
		{	const auto arg = 1 + (vnew - vold) / vt;
			if (arg > 0)
				return vold + vt * log(arg);
			else
				return vcrit;
		}
		else
			return vt *log(vnew/vt);
	}
	else
		return vnew;
#endif
#endif
	}
#endif
	auto calculate(
		const std::array<double, std::size_t(enumCircuitNodes::NumberOfNodes)>& _rNodeVoltages,
		lufac::index2Index2Index2Double&_r2,
		lufac::index2Index2Double&_rO,
		lufac::index2Double &_rV,
		const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> & _pT
	) const
	{
#define __create__(Vbe, b, e) const auto Vbe = createIndep<enumBranches::Vbe>(\
	enumNodes::b != enumNodes::NumberOfNodes && _pT[std::size_t(enumNodes::b)] != std::size_t(enumCircuitNodes::NumberOfNodes)\
		? (enumNodes::e != enumNodes::NumberOfNodes && _pT[std::size_t(enumNodes::e)] != std::size_t(enumCircuitNodes::NumberOfNodes)\
			? _rNodeVoltages[_pT[std::size_t(enumNodes::b)]] - _rNodeVoltages[_pT[std::size_t(enumNodes::e)]]\
			: _rNodeVoltages[_pT[std::size_t(enumNodes::b)]]\
		)\
		: (enumNodes::e != enumNodes::NumberOfNodes && _pT[std::size_t(enumNodes::e)] != std::size_t(enumCircuitNodes::NumberOfNodes)\
			? - _rNodeVoltages[_pT[std::size_t(enumNodes::e)]]\
			: 0.0\
		),\
		false\
	);
#define __COMMA__
#include "inputs.h"
#ifdef __DIODE__
		static constexpr auto IS = 1e-12;
		constexpr auto VT = 0.025;
		constexpr auto R1 = 1e3;
#if 0
		const auto ID0 = IS*(exp(VD0/VT) - 1.0);
#else
		const auto ID0 = IS*(EXP(VD0/VT, _VD0/VT, _pbLimit) - 1.0);
#endif
		const auto IR0 = VR0/R1;
#if 0
		const auto ID1 = IS*(exp(VD1/VT) - 1.0);
#else
		const auto ID1 = IS*(EXP(VD1/VT, _VD1/VT, _pbLimit) - 1.0);
#endif
		const auto IR1 = VR1/R1;
#else
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
	const auto Itrp  = ISP * ( exp ( Vbcp / ( NFP * Vtv ) ) - 1.0 );
	const auto q2p   = Itfp * IIKP;				// only fwd part
	const auto qbp   = 0.5 * ( 1.0 + sqrt ( 1.0 + 4.0 * q2p ) );// no Early effect
	const auto Iccp  = ( Itfp - Itrp ) / qbp;			// parasitic transport I

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
			assert(std::isfinite(value(Vbcx)) && !std::isnan(value(Vbcx)));
			const auto Kbci = sqrt ( 1.0 + GAMM * exp ( Vbci / Vtv ) );
			assert(std::isfinite(value(Kbci)) && !std::isnan(value(Kbci)));
			const auto Kbcx = sqrt ( 1.0 + GAMM * exp ( Vbcx / Vtv ) );
			assert(std::isfinite(value(Kbcx)) && !std::isnan(value(Kbcx)));
			const auto rKp1 = ( Kbci + 1.0 ) / ( Kbcx + 1.0 );
			assert(std::isfinite(value(rKp1)) && !std::isnan(value(rKp1)));
			const auto Iohm = ( Vrci + Vtv * ( Kbci - Kbcx - log ( rKp1 ) ) ) / RCI;
			assert(std::isfinite(value(Iohm)) && !std::isnan(value(Iohm)));
			const auto derf = IVO * RCI * Iohm / ( 1.0 + 0.5 * IVO * IHRCF *
				sqrt ( Vrci * Vrci +  0.01 ) );
			assert(std::isfinite(value(derf)) && !std::isnan(value(derf)));
			const auto Irci = Iohm / sqrt ( 1.0 + derf * derf  );
			assert(std::isfinite(value(Irci)) && !std::isnan(value(Irci)));
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
		{
			const auto Ith   =   Ibe  * Vbei + Ibc  * Vbci
			+ ( Itzf - Itzr ) * ( Vbei - Vbci )
			+ Ibep * Vbep + Ibcp * Vbcp + Iccp * ( Vbep - Vbcp )
			+ Ircx * Vrcx + Irci * Vrci + Irbx * Vrbx
			+ Irbi * Vrbi + Ire  * Vre  + Irbp * Vrbp
			+ Irs  * Vrs  + Ibex * Vbex - Igc  * Vbci;

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
#endif
#define __create__(a, b, c)\
{	{	const auto iBT = _pT[std::size_t(enumNodes::b)];\
		if (iBT != std::size_t(enumCircuitNodes::NumberOfNodes))\
			writeOutput(_rV[iBT], _rO[iBT], _r2[iBT], -a, _pT);\
	}\
	{	const auto iCT = _pT[std::size_t(enumNodes::c)];\
		if (iCT != std::size_t(enumCircuitNodes::NumberOfNodes))\
			writeOutput(_rV[iCT], _rO[iCT], _r2[iCT], a, _pT);\
	}\
}
#define __COMMA__
#include "currentSources.h"
//	End of equations.
	}
	template<typename T>
	struct copyOutput
	{	lufac::index2Double&m_rO;
		lufac::index2Index2Double&m_r1;
		double&m_rValue;
		const ctaylor<T, 2>&m_rV;
		const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> &m_pT;
		copyOutput(
			double&_rValue,
			lufac::index2Double&_rO,
			lufac::index2Index2Double&_r1,
			const ctaylor<T, 2>&_rV,
			const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> &_pT
		)
			:m_rO(_rO),
			m_r1(_r1),
			m_rV(_rV),
			m_rValue(_rValue),
			m_pT(_pT)
		{
		}
		constexpr enumCircuitNodes translate(const enumNodes _e) const
		{	if (_e == enumNodes::NumberOfNodes)
				return enumCircuitNodes::NumberOfNodes;
			else
				return enumCircuitNodes(m_pT[std::size_t(_e)]);
		}
		constexpr std::pair<enumCircuitNodes, enumCircuitNodes> translate(const std::pair<enumNodes, enumNodes>&_r) const
		{	return std::make_pair(translate(_r.first), translate(_r.second));
		}
		template<std::size_t POS>
		void operator()(const enumCircuitNodes _e0, const enumCircuitNodes _e1, const mp_size_t<POS>&, const bool _bMinus, const bool _bDouble = false) const
		{	if (_e0 != enumCircuitNodes::NumberOfNodes && _e1 != enumCircuitNodes::NumberOfNodes)
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
		template<typename ENUM>
		void operator()(const mp_list<mp_list<ENUM, mp_size_t<2> > >&) const
		{	const auto sNodePair = translate(s_aInput2NodePair[ENUM::value]);
			typedef mp_list<mp_list<ENUM, mp_size_t<1> > > LIST;
			constexpr std::size_t POS = mp_find<T, LIST>::value;
			(*this)(sNodePair.first, sNodePair.first, mp_size_t<POS>(), false, true);
			(*this)(sNodePair.second, sNodePair.second, mp_size_t<POS>(), false, true);
			(*this)(sNodePair.first, sNodePair.second, mp_size_t<POS>(), true, true);
			(*this)(sNodePair.second, sNodePair.first, mp_size_t<POS>(), true, true);
		}
		template<typename ENUM0, typename ENUM1>
		void operator()(const mp_list<mp_list<ENUM0, mp_size_t<1> >, mp_list<ENUM1, mp_size_t<1> > >&) const
		{	const auto sNodePair0 = translate(s_aInput2NodePair[ENUM0::value]);
			const auto sNodePair1 = translate(s_aInput2NodePair[ENUM1::value]);
			typedef mp_list<mp_list<ENUM0, mp_size_t<1> >, mp_list<ENUM1, mp_size_t<1> > > LIST;
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
		template<typename ENUM>
		void operator()(const mp_list<mp_list<ENUM, mp_size_t<1> > >&) const
		{	const auto sNodePair = translate(s_aInput2NodePair[ENUM::value]);
			typedef mp_list<mp_list<ENUM, mp_size_t<1> > > LIST;
			constexpr std::size_t POS = mp_find<T, LIST>::value;
			if (const double d = m_rV.m_s.at(POS))
			{	if (sNodePair.first != enumCircuitNodes::NumberOfNodes)
					m_rO[std::size_t(sNodePair.first)] -= d;
				if (sNodePair.second != enumCircuitNodes::NumberOfNodes)
					m_rO[std::size_t(sNodePair.second)] += d;
			}
		}
		void operator()(const mp_list<>&) const
		{	m_rValue += m_rV.m_s.at(mp_find<T, mp_list<> >::value);
		}
	};
	template<typename T>
	static void writeOutput(double&_rValue, lufac::index2Double&_rO, lufac::index2Index2Double&_r1, const ctaylor<T, 2>&_rV, const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> & _pT)
	{	mp_for_each<T>(copyOutput<T>(_rValue, _rO, _r1, _rV, _pT));
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
constexpr enumCircuitNodes translateNodes(const enumNodes _e)
{	switch (_e)
	{	default:
			throw std::logic_error("Invalid node ID!");
		case enumNodes::NumberOfNodes:
			return enumCircuitNodes::NumberOfNodes;
#ifdef __DIODE__
		case enumNodes::a0:
			return enumCircuitNodes::a0;
		case enumNodes::b0:
			return enumCircuitNodes::b0;
		case enumNodes::c0:
			return enumCircuitNodes::NumberOfNodes;
		case enumNodes::a1:
			return enumCircuitNodes::a1;
		case enumNodes::b1:
			return enumCircuitNodes::b1;
		case enumNodes::c1:
			return enumCircuitNodes::NumberOfNodes;
#else
		case enumNodes::c:
			return enumCircuitNodes::c;
		case enumNodes::b:
			return enumCircuitNodes::b;
		case enumNodes::e:
			return enumCircuitNodes::e;
		case enumNodes::s:
			return enumCircuitNodes::s;
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
#endif
	}
}
lufac::index2Index2Double transpose(const lufac::index2Index2Double&_r)
{	lufac::index2Index2Double s;
	for (const auto &rR : _r)
		for (const auto &rC : rR.second)
			s[rC.first][rR.first] = rC.second;
	return s;
}
lufac::index2Index2Double sqrT2(const lufac::index2Index2Double&_r)
{	lufac::index2Index2Double sRet;
	const auto sT = transpose(_r);
	for (const auto &rR0 : _r)
	{	auto &rS = sRet[rR0.first];
		for (const auto &rC1 : sT)
		{	auto &rT = rS[rC1.first];
			for (auto p0 = rR0.second.begin(), p1 = rC1.second.begin(); p0 != rR0.second.end() && p1 != rC1.second.end();)
				if (p0->first < p1->first)
					++p0;
				else
				if (p0->first > p1->first)
					++p1;
				else
					rT += p0++->second*p1++->second;
			if (rT == 0.0)
				rS.erase(rC1.first);
			else
				rT *= 2.0;
		}
	}
	return sRet;
}
lufac::index2Index2Double FtF2(const lufac::index2Double&_rF, const lufac::index2Index2Index2Double&_rF2)
{	lufac::index2Index2Double sRet;
	for (const auto &r0 : _rF2)
	{	auto &rT0 = sRet[r0.first];
		for (const auto &r1 : r0.second)
		{	auto &rT1 = rT0[r1.first];
			for (const auto &r2 : r1.second)
				rT1 += r2.second*_rF.at(r2.first);
		}
	}
	return sRet;
}
lufac::index2Double twoFtF1(const lufac::index2Double&_rF, const lufac::index2Index2Double&_rF1)
{	lufac::index2Double sRet;
	for (const auto &rR : _rF1)
	{	auto &rT = sRet[rR.first];
		for (const auto &rC : rR.second)
			rT += rC.second*_rF.at(rC.first);
		rT *= 2.0;
	}	
	return sRet;
}
lufac::index2Double operator-(const lufac::index2Double&_r0, const lufac::index2Double&_r1)
{	lufac::index2Double s;
	for (auto pR0 = _r0.begin(), pR1 = _r1.begin(); pR0 != _r0.end() && pR1 != _r1.end(); )
		if (pR0->first < pR1->first)
		{	s[pR0->first] = pR0->second;
			++pR0;
		}
		else
		if (pR0->first > pR1->first)
		{	s[pR1->first] = -pR1->second;
			++pR1;
		}
		else
		{	s[pR1->first] = pR0->second - pR1->second;
			++pR0;
			++pR1;
		}
	return s;
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
lufac::index2Double operator-(const lufac::index2Double&_r)
{	lufac::index2Double s;
	for (const auto &r : _r)
		s[r.first] = -r.second;
	return s;
}
#if 0
lufac::index2Index2Double operator-(const lufac::index2Index2Double&_r0, const lufac::index2Index2Double&_r1)
{	lufac::index2Index2Double sRet;
	lufac::index2Index2Double::const_iterator pR0, pR1;
	for (pR0 = _r0.begin(), pR1 = _r1.begin(); pR0 != _r0.end() && pR1 != _r1.end(); )
		if (pR0->first < pR1->first)
		{	sRet[pR0->first] = pR0->second;
			++pR0;
		}
		else
		if (pR0->first > pR1->first)
		{	sRet[pR1->first] = -pR1->second;
			++pR1;
		}
		else
		{	sRet[pR0->first] = pR0->second - pR1->second;
			++pR0;
			++pR1;
		}
	for (; pR0 != _r0.end(); ++pR0)
		sRet[pR0->first] = pR0->second;
	for (; pR1 != _r1.end(); ++pR1)
		sRet[pR1->first] = -pR1->second;
	return sRet;
}
#endif
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
}
int main(int argc, char**argv)
try
{	using namespace vbic95;
	using namespace lufac;
	if (argc < 2)
	{	std::cerr << argv[0] << ": usage: " << argv[0] << " PARS" << std::endl;
		return 1;
	}
	vbic sI(readParams(argv[1]));
	const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> sTrans{
#define __create__(a) std::size_t(translateNodes(enumNodes::a))
#define __COMMA__ ,
#include "nodes.h"
	};
#ifndef __DIODE__
	for (double vb = 0.7; vb <= 1.00001; vb += 0.05)
		for (double vc = 0.0; vc <= 5.00001; vc += 0.05)
#else
		for (double VD = 0.0; VD <= 5.00001; VD += 0.05)
#endif
		{	std::array<double, std::size_t(enumCircuitNodes::NumberOfNodes)> sV({});
			std::array<double, std::size_t(enumCircuitNodes::NumberOfNodes)> sV1({});
#ifdef __DIODE__
			sV[std::size_t(enumCircuitNodes::b0)] = VD;
			sV[std::size_t(enumCircuitNodes::b1)] = VD;
#else
			sV[std::size_t(enumCircuitNodes::b)] = vb;
			sV[std::size_t(enumCircuitNodes::bx)] = vb;
			sV[std::size_t(enumCircuitNodes::bi)] = vb;
			sV[std::size_t(enumCircuitNodes::bp)] = vc;
			sV[std::size_t(enumCircuitNodes::e)] = 0.0;
			sV[std::size_t(enumCircuitNodes::ei)] = 0.0;
			sV[std::size_t(enumCircuitNodes::s)] = 0.0;
			sV[std::size_t(enumCircuitNodes::si)] = 0.0;
			sV[std::size_t(enumCircuitNodes::c)] = vc;
			sV[std::size_t(enumCircuitNodes::ci)] = vc;
			sV[std::size_t(enumCircuitNodes::cx)] = vc;
#endif
	//ve=0;vs=0
			//for(vb=0.7;vb<=1.00001;vb+=0.05)
			//for(vc=0.0;vc<=5.00001;vc+=0.05)
		//print vc,vb,ve,vs
			for (std::size_t i = 0; true; ++i)
			{	index2Index2Index2Double sH;
				index2Index2Double sJ;
				index2Double sValues;
				sI.calculate(sV, sH, sJ, sValues, sTrans);
				for (const auto &rS : sJ)
					for (const auto &r : rS.second)
						if (!std::isfinite(r.second) || std::isnan(r.second))
							throw std::logic_error("Calculate produces NAN!");
				for (const auto &r : sValues)
					if (!std::isfinite(r.second) || std::isnan(r.second))
						throw std::logic_error("Calculate produces NAN!");
#ifdef __DIODE__
				sValues[std::size_t(enumCircuitNodes::IVB0)] += VD - sV.at(std::size_t(enumCircuitNodes::b0));
				sValues[std::size_t(enumCircuitNodes::b0)] += - sV.at(std::size_t(enumCircuitNodes::IVB0));
				sJ[std::size_t(enumCircuitNodes::b0)][std::size_t(enumCircuitNodes::IVB0)] += 1.0;
				sJ[std::size_t(enumCircuitNodes::IVB0)][std::size_t(enumCircuitNodes::b0)] += 1.0;
				
				sValues[std::size_t(enumCircuitNodes::IVB1)] += VD - sV.at(std::size_t(enumCircuitNodes::b1));
				sValues[std::size_t(enumCircuitNodes::b1)] += - sV.at(std::size_t(enumCircuitNodes::IVB1));
				sJ[std::size_t(enumCircuitNodes::b1)][std::size_t(enumCircuitNodes::IVB1)] += 1.0;
				sJ[std::size_t(enumCircuitNodes::IVB1)][std::size_t(enumCircuitNodes::b1)] += 1.0;
#else
				sValues[std::size_t(enumCircuitNodes::IVC)] += vc - sV.at(std::size_t(enumCircuitNodes::c));
				sValues[std::size_t(enumCircuitNodes::IVB)] += vb - sV.at(std::size_t(enumCircuitNodes::b));
				sValues[std::size_t(enumCircuitNodes::IVE)] += - sV.at(std::size_t(enumCircuitNodes::e));
				sValues[std::size_t(enumCircuitNodes::IVS)] += - sV.at(std::size_t(enumCircuitNodes::s));
				
				sValues[std::size_t(enumCircuitNodes::c)] += -sV.at(std::size_t(enumCircuitNodes::IVC));
				sValues[std::size_t(enumCircuitNodes::b)] += -sV.at(std::size_t(enumCircuitNodes::IVB));
				sValues[std::size_t(enumCircuitNodes::e)] += -sV.at(std::size_t(enumCircuitNodes::IVE));
				sValues[std::size_t(enumCircuitNodes::s)] += -sV.at(std::size_t(enumCircuitNodes::IVS));
				
				sJ[std::size_t(enumCircuitNodes::c)][std::size_t(enumCircuitNodes::IVC)] += 1.0;
				sJ[std::size_t(enumCircuitNodes::b)][std::size_t(enumCircuitNodes::IVB)] += 1.0;
				sJ[std::size_t(enumCircuitNodes::e)][std::size_t(enumCircuitNodes::IVE)] += 1.0;
				sJ[std::size_t(enumCircuitNodes::s)][std::size_t(enumCircuitNodes::IVS)] += 1.0;
				
				sJ[std::size_t(enumCircuitNodes::IVC)][std::size_t(enumCircuitNodes::c)] += 1.0;
				sJ[std::size_t(enumCircuitNodes::IVB)][std::size_t(enumCircuitNodes::b)] += 1.0;
				sJ[std::size_t(enumCircuitNodes::IVE)][std::size_t(enumCircuitNodes::e)] += 1.0;
				sJ[std::size_t(enumCircuitNodes::IVS)][std::size_t(enumCircuitNodes::s)] += 1.0;
#endif
				const auto sM = sqrT2(sJ) + FtF2(sValues, sH);
				auto sValues1 = twoFtF1(sValues, sJ);
				const auto sFactored = factor(sM, sValues1);
				const auto sDelta = solve(sFactored, sValues1);
				for (const auto &r : sDelta)
					if (!std::isfinite(r.second) || std::isnan(r.second))
						throw std::logic_error("solve produces NAN!");
				const auto dNormO = std::sqrt(
					std::accumulate(
						sValues1.begin(),
						sValues1.end(),
						0.0,
						[](const double _dSum, const index2Double::value_type&_r)
						{	return _dSum + _r.second*_r.second;
						}
					)/sValues1.size()
				);
				const auto dNorm = std::sqrt(
					std::accumulate(
						sDelta.begin(),
						sDelta.end(),
						0.0,
						[](const double _dSum, const index2Double::value_type&_r)
						{	return _dSum + _r.second*_r.second;
						}
					)/sDelta.size()
				);
				for (const auto &r : sDelta)
					sV[r.first] += r.second;
				if (dNorm < 1e-9 && dNormO < 1e-9)
				{	std::cout << "iter=" << i + 1 << ',';
					break;
				}
			}
#ifdef __DIODE__
			for (const auto i : {enumCircuitNodes::b0, enumCircuitNodes::IVB0, enumCircuitNodes::IVB1})
#else
			for (const auto i : {enumCircuitNodes::c, enumCircuitNodes::b, enumCircuitNodes::IVC, enumCircuitNodes::IVB, enumCircuitNodes::IVS})
#endif
				std::cout << sV[std::size_t(i)] << ",";
			std::cout << "\n";
		}
} catch (const std::exception&_r)
{	std::cerr << argv[1] << ": Exception caught: " << _r.what() << std::endl;
	return 1;
}
