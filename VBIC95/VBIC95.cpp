#include "../ctaylor.h"
#include "../LUFAC/lufac.h"
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <numeric>
//#define SELF_HEATING
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
static enumMembers2double readParams(void)
{	std::ifstream sFile("PARS");
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

enum class eCircuitNodes
{	c,
	b,
	e,
	s,
#ifdef SELF_HEATING
	dt,
//	tl,
#endif
	cx,
	ci,
	bx,
	bi,
	ei,
	si,
	bp,
#ifdef EXCESS_PHASE
	xf1,
	xf2,
#endif
	IVC,
	IVB,
	IVE,
	IVS,
	NumberOfNodes
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
#include "outputs.h"
};
template<enumBranches E>
using createIndep = ctaylor<makeIndependent<std::size_t(E)>, 1>;
struct vbic
{
#define __create__(a) const double a;
#define __COMMA__
#include "members.h"
	vbic(const enumMembers2double &_r = readParams())
		:
#define __create__(a) a(_r.at(enumMembers::a))
#define __COMMA__ ,
#include "members.h"
	{
	}
	auto calculate(
		const std::array<double, std::size_t(eCircuitNodes::NumberOfNodes)>& _rNodeVoltages,
		lufac::index2Index2Double&_rO,
		lufac::index2Double &_rV,
		const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> & _pT
	) const
	{
#define __create__(Vbe, b, e) const auto Vbe = createIndep<enumBranches::Vbe>(\
	enumNodes::b != enumNodes::NumberOfNodes && _pT[std::size_t(enumNodes::b)] != std::size_t(eCircuitNodes::NumberOfNodes)\
		? (enumNodes::e != enumNodes::NumberOfNodes && _pT[std::size_t(enumNodes::e)] != std::size_t(eCircuitNodes::NumberOfNodes)\
			? _rNodeVoltages[_pT[std::size_t(enumNodes::b)]] - _rNodeVoltages[_pT[std::size_t(enumNodes::e)]]\
			: _rNodeVoltages[_pT[std::size_t(enumNodes::b)]]\
		)\
		: (enumNodes::e != enumNodes::NumberOfNodes && _pT[std::size_t(enumNodes::e)] != std::size_t(eCircuitNodes::NumberOfNodes)\
			? - _rNodeVoltages[_pT[std::size_t(enumNodes::e)]]\
			: 0.0\
		),\
		false\
	);
#define __COMMA__
#include "inputs.h"
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
			const auto Kbci = sqrt ( 1.0 + GAMM * exp ( Vbci / Vtv ) );
			const auto Kbcx = sqrt ( 1.0 + GAMM * exp ( Vbcx / Vtv ) );
			const auto rKp1 = ( Kbci + 1.0 ) / ( Kbcx + 1.0 );
			const auto Iohm = ( Vrci + Vtv * ( Kbci - Kbcx - log ( rKp1 ) ) ) / RCI;
			const auto derf = IVO * RCI * Iohm / ( 1.0 + 0.5 * IVO * IHRCF *
				sqrt ( Vrci * Vrci +  0.01 ) );
			const auto Irci = Iohm / sqrt ( 1.0 + derf * derf  );
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

#define __create__(a, b, c)\
{	if (enumNodes::b != enumNodes::NumberOfNodes)\
	{	const auto iBT = _pT[std::size_t(enumNodes::b)];\
		if (iBT != std::size_t(eCircuitNodes::NumberOfNodes))\
			writeOutput(_rV[iBT], _rO[iBT], -a, _pT);\
	}\
	if (enumNodes::c != enumNodes::NumberOfNodes)\
	{	const auto iCT = _pT[std::size_t(enumNodes::c)];\
		if (iCT != std::size_t(eCircuitNodes::NumberOfNodes))\
			writeOutput(_rV[iCT], _rO[iCT], a, _pT);\
	}\
}
#define __COMMA__
#include "outputs.h"
//	End of equations.
	}
	template<typename T>
	struct copyOutput
	{	lufac::index2Double&m_rO;
		double&m_rValue;
		const ctaylor<T, 1>&m_rV;
		const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> &m_pT;
		copyOutput(
			double&_rValue,
			lufac::index2Double&_rO,
			const ctaylor<T, 1>&_rV,
			const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> &_pT
		)
			:m_rO(_rO),
			m_rV(_rV),
			m_rValue(_rValue),
			m_pT(_pT)
		{
		}
		template<typename ENUM>
		void operator()(const mp_list<mp_list<ENUM, mp_size_t<1> > >&) const
		{	const auto &rNodePair = s_aInput2NodePair[ENUM::value];
			typedef mp_list<mp_list<ENUM, mp_size_t<1> > > LIST;
			constexpr std::size_t POS = mp_find<T, LIST>::value;
			if (rNodePair.first != enumNodes::NumberOfNodes)
			{	const auto iT = m_pT[std::size_t(rNodePair.first)];
				if (iT != std::size_t(eCircuitNodes::NumberOfNodes))
					m_rO[iT] -= m_rV.m_s.at(POS);
			}
			if (rNodePair.second != enumNodes::NumberOfNodes)
			{	const auto iT = m_pT[std::size_t(rNodePair.second)];
				if (iT != std::size_t(eCircuitNodes::NumberOfNodes))
					m_rO[iT] += m_rV.m_s.at(POS);
			}
		}
		void operator()(const mp_list<>&) const
		{	constexpr std::size_t POS = mp_find<T, mp_list<> >::value;
			m_rValue += m_rV.m_s.at(POS);
		}
		template<typename LIST>
		void operator()(const LIST&) const
		{
		}
	};
	template<typename T>
	static void writeOutput(double&_rValue, lufac::index2Double&_rO, const ctaylor<T, 1>&_rV, const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> & _pT)
	{	mp_for_each<T>(copyOutput<T>(_rValue, _rO, _rV, _pT));
	}
	static const char*const s_aNames[];
};
const char*const vbic::s_aNames[] = {
#define __create__(a, b, c) #a
#define __COMMA__ ,
#include "outputs.h"
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
constexpr eCircuitNodes translateNodes(const enumNodes _e)
{	switch (_e)
	{	default:
			throw std::logic_error("Invalid node ID!");
		case enumNodes::NumberOfNodes:
			return eCircuitNodes::NumberOfNodes;
		case enumNodes::c:
			return eCircuitNodes::c;
		case enumNodes::b:
			return eCircuitNodes::b;
		case enumNodes::e:
			return eCircuitNodes::e;
		case enumNodes::s:
			return eCircuitNodes::s;
#ifdef SELF_HEATING
		case enumNodes::dt:
			return eCircuitNodes::dt;
		case enumNodes::tl:
			return eCircuitNodes::NumberOfNodes;
#endif
		case enumNodes::cx:
			return eCircuitNodes::cx;
		case enumNodes::ci:
			return eCircuitNodes::ci;
		case enumNodes::bx:
			return eCircuitNodes::bx;
		case enumNodes::bi:
			return eCircuitNodes::bi;
		case enumNodes::ei:
			return eCircuitNodes::ei;
		case enumNodes::si:
			return eCircuitNodes::si;
		case enumNodes::bp:
			return eCircuitNodes::bp;
#ifdef EXCESS_PHASE
		case enumNodes::xf1:
			return eCircuitNodes::xf1;
		case enumNodes::xf2:
			return eCircuitNodes::xf2;
#endif
	}
}
}
int main(int, char**)
{	using namespace vbic95;
	using namespace lufac;
	vbic sI;
	const std::array<std::size_t, std::size_t(1) + std::size_t(enumNodes::NumberOfNodes)> sTrans{
#define __create__(a) std::size_t(translateNodes(enumNodes::a))
#define __COMMA__ ,
#include "nodes.h"
	};
	for (double vb = 0.7; vb <= 1.00001; vb += 0.05)
		for (double vc = 0.0; vc <= 5.00001; vc += 0.05)
		{	std::array<double, std::size_t(eCircuitNodes::NumberOfNodes)> sV({});
			sV[std::size_t(eCircuitNodes::b)] = vb;
			sV[std::size_t(eCircuitNodes::bx)] = vb;
			sV[std::size_t(eCircuitNodes::bi)] = vb;
			sV[std::size_t(eCircuitNodes::bp)] = vc;
			sV[std::size_t(eCircuitNodes::e)] = 0.0;
			sV[std::size_t(eCircuitNodes::ei)] = 0.0;
			sV[std::size_t(eCircuitNodes::s)] = 0.0;
			sV[std::size_t(eCircuitNodes::si)] = 0.0;
			sV[std::size_t(eCircuitNodes::c)] = vc;
			sV[std::size_t(eCircuitNodes::ci)] = vc;
			sV[std::size_t(eCircuitNodes::cx)] = vc;

	//{ve=0;vs=0
		 //for(vb=0.7;vb<=1.00001;vb+=0.05){
		 //for(vc=0.0;vc<=5.00001;vc+=0.05){
		//print vc,vb,ve,vs
			for (std::size_t i = 0; i < 100; ++i)
			{	index2Index2Double s;
				index2Double sValues;
				sI.calculate(sV, s, sValues, sTrans);
				
				sValues[std::size_t(eCircuitNodes::IVC)] += vc - sV.at(std::size_t(eCircuitNodes::c));
				sValues[std::size_t(eCircuitNodes::IVB)] += vb - sV.at(std::size_t(eCircuitNodes::b));
				sValues[std::size_t(eCircuitNodes::IVE)] += - sV.at(std::size_t(eCircuitNodes::e));
				sValues[std::size_t(eCircuitNodes::IVS)] += - sV.at(std::size_t(eCircuitNodes::s));
				
				sValues[std::size_t(eCircuitNodes::c)] += -sV.at(std::size_t(eCircuitNodes::IVC));
				sValues[std::size_t(eCircuitNodes::b)] += -sV.at(std::size_t(eCircuitNodes::IVB));
				sValues[std::size_t(eCircuitNodes::e)] += -sV.at(std::size_t(eCircuitNodes::IVE));
				sValues[std::size_t(eCircuitNodes::s)] += -sV.at(std::size_t(eCircuitNodes::IVS));
				
				s[std::size_t(eCircuitNodes::c)][std::size_t(eCircuitNodes::IVC)] += 1.0;
				s[std::size_t(eCircuitNodes::b)][std::size_t(eCircuitNodes::IVB)] += 1.0;
				s[std::size_t(eCircuitNodes::e)][std::size_t(eCircuitNodes::IVE)] += 1.0;
				s[std::size_t(eCircuitNodes::s)][std::size_t(eCircuitNodes::IVS)] += 1.0;
				
				s[std::size_t(eCircuitNodes::IVC)][std::size_t(eCircuitNodes::c)] += 1.0;
				s[std::size_t(eCircuitNodes::IVB)][std::size_t(eCircuitNodes::b)] += 1.0;
				s[std::size_t(eCircuitNodes::IVE)][std::size_t(eCircuitNodes::e)] += 1.0;
				s[std::size_t(eCircuitNodes::IVS)][std::size_t(eCircuitNodes::s)] += 1.0;
				const auto sFactored = factor(s, sValues);
				const auto sDelta = solve(sFactored, sValues);
				const auto dNormO = std::sqrt(
					std::accumulate(
						sValues.begin(),
						sValues.end(),
						0.0,
						[](const double _dSum, const index2Double::value_type&_r)
						{	return _dSum + _r.second*_r.second;
						}
					)/sValues.size()
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
					break;
			}
			for (const auto i : {eCircuitNodes::c, eCircuitNodes::b, eCircuitNodes::IVC, eCircuitNodes::IVB, eCircuitNodes::IVS})
				std::cout << sV[std::size_t(i)] << ",";
			std::cout << "\n";
		}
}
