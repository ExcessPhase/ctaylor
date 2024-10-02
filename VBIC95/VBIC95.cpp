#include "../ctaylor.h"
#include <map>
#include <string>
#include <fstream>
//#include <initializer_list>
#define SELF_HEATING
namespace vbic95
{
typedef std::map<std::string, double> string2double;
static string2double readParams(void)
{	std::ifstream sFile("PARS");
	double d;
	std::string s;
	string2double sRet;
	while (true)
	{	if (sFile >> d >> s)
		{	sRet[s] = d;
			sRet[s + "_TNOM"] = d;
		}
		else
			break;
	}
	return sRet;
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

enum class enumNodes
{	c,
	b,
	e,
	s,
#ifdef SELF_HEATING
	dt,
	tl,
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
	NumberOfNodes,
};
enum class enumBranches:std::size_t
{	Vbe,
	Vbc,
	Vbei,
	Vbex,
	Vbci,
	Vrcx,
	Vrci,
	Vrbx,
	Vrbi,
	Vre,
	Vrs,
	Vbep,
	Vbcp,
	Vrbp,
#ifdef SELF_HEATING
	delT,
#endif
#ifdef EXCESS_PHASE
	Vcxf,
	Vrxf,
#endif
};
enum class enumCurrentOutputs:std::size_t
{
#define __create__(a, b, c) a
#define __COMMA__ ,
#include "outputs.h"
};
template<enumBranches E>
using createIndep = ctaylor<makeIndependent<std::size_t(E)>, 2>;
static constexpr double KB = 1.380662E-23;	// Boltzmann's constant, J/K
static constexpr double QQ = 1.602189E-19;	// magnitude of electronic charge, C
static constexpr double TABS = 2.731500E+02;	// 0 degrees C in degrees K
static constexpr double TAMB = 27.0;
struct vbic
{
#define __create__(a) const double a;
#define __COMMA__
#include "members.h"
#if 0
	vbic(void)
		:
#define __create__(a) a(),
#include "members.h"
		COMMA()
	{
	}
#endif
	vbic(const string2double&_r = readParams())
		:
#define __create__(a) a(_r.at(#a))
#define __COMMA__ ,
#include "members.h"
	{
	}
	auto calculate(const std::array<double, std::size_t(enumNodes::NumberOfNodes)>& _rNodeVoltages) const
	{
#define createVoltage(Vbe, b, e) const auto Vbe = createIndep<enumBranches::Vbe>(_rNodeVoltages[std::size_t(enumNodes::b)] - _rNodeVoltages[std::size_t(enumNodes::e)], false)
		createVoltage(Vbe, b, e);
		createVoltage(Vbc, b, c);
		createVoltage(Vbei, bi, ei);	// intrinsic b-e voltage
		createVoltage(Vbex, bx, ei);	// side      b-e voltage
		createVoltage(Vbci, bi, ci);	// intrinsic b-c voltage
		createVoltage(Vrcx, c, cx);	// voltage across RCX
		createVoltage(Vrci, cx, ci);	// voltage across RCI
		createVoltage(Vrbx, b, bx);	// voltage across RBX
		createVoltage(Vrbi, bx, bi);	// voltage across RBI
		createVoltage(Vre, e, ei);	// voltage across RE
		createVoltage(Vrs, s, si);	// voltage across RS
		createVoltage(Vbep, bx, bp);	// parasitic b-e voltage (pnp polarity)
		createVoltage(Vbcp, si, bp);	// parasitic b-c voltage (pnp polarity)
		createVoltage(Vrbp, cx, bp);	// voltage across RBP
#ifdef SELF_HEATING
		createVoltage(delT, dt, tl);// voltage across RTH, local temperature rise
					// measured with respect to ground (ambient)
#endif
#ifdef EXCESS_PHASE
#define createVoltage1(Vbe, b) const auto Vbe = createIndep<enumBranches::Vbe>(_rNodeVoltages[enumNodes::b], false)
		createVoltage1(Vcxf, xf1);	// voltage across excess-phase capacitor
		createVoltage1(Vrxf, xf2);	// voltage across excess-phase resistor Itxf
#endif
#ifdef SELF_HEATING
//	This section defines mappings of temperature dependent parameters.
//	Note that with SELF-HEATING these mappings must be done at each
//	bias calculation, and cannot be relegated to preprocessing, because
//	the "voltage" for the local temperature rise delT must be taken
//	into account.

#ifdef SELF_HEATING
	const auto Tdev  = TAMB + TABS + delT;	// device temperature (K), with delT
#else
	const auto Tdev  = TAMB + TABS;		// device temperature (K)
#endif
	const auto Tini  = TNOM + TABS;		// TNOM in K

	const auto Vtv   = KB * Tdev / QQ;		// thermal voltage
	const auto rT    = Tdev / Tini;		// ratio Tdev/TNOM in K

	const auto RCX   = RCX_TNOM * pow(rT, XRC);	// temperature mapping
	const auto RCI   = RCI_TNOM * pow(rT, XRC);	// of resistances is
	const auto RBX   = RBX_TNOM * pow(rT, XRB);	// according to mobility
	const auto RBI   = RBI_TNOM * pow(rT, XRB);	// variation, which
	const auto RE    = RE_TNOM  * pow(rT, XRE);	// should be same in
	const auto RS    = RS_TNOM  * pow(rT, XRS);	// emit/coll & base/subs
	const auto RBP   = RBP_TNOM * pow(rT, XRC);	// Note: RBP is in coll

//	Note: the following differs from the standard SPICE temperature
//	mappings for IS/ISE/BF. This is for two reasons. First, 
//	the base current is formulated directly in terms
//	of ideal and non-ideal currents rather than in terms of current
//	gain BF/transport current and non-ideal current. Second,
//	particularly for HBTs, using the bandgap in the first-order
//	theory expression for how IS should change with temperature often
//	gives a poor model. To properly track IS and beta (both low and
//	moderate bias) over temperature there is an activation energy
//	in the first-order model that is (slightly) different from the bandgap,
//	and is different for all of IS, IBEI and IBEN.
//	01/26/94: the b-e and b-c components show slightly different behavior
//	over temperature, so "activation energies" have been introduced for
//	ideal and non-ideal components of b-e, b-c and s-c junctions. This
//	allows separate fitting of forward and reverse beta curves over
//	temperature.

	const auto IS    = IS_TNOM    * pow( pow(rT, XIS) * exp ( - EA   * ( 1.0 - rT ) / Vtv ), 1.0 / NF_TNOM );
	const auto ISP   = ISP_TNOM   * pow( pow(rT, XIS) * exp ( - EA   * ( 1.0 - rT ) / Vtv ), 1.0 / NFP  );
	const auto IBEI  = IBEI_TNOM  * pow( pow(rT, XII) * exp ( - EAIE * ( 1.0 - rT ) / Vtv ), 1.0 / NEI  );
	const auto IBEN  = IBEN_TNOM  * pow( pow(rT, XIN) * exp ( - EANE * ( 1.0 - rT ) / Vtv ), 1.0 / NEN  );
	const auto IBCI  = IBCI_TNOM  * pow( pow(rT, XII) * exp ( - EAIC * ( 1.0 - rT ) / Vtv ), 1.0 / NCI  );
	const auto IBCN  = IBCN_TNOM  * pow( pow(rT, XIN) * exp ( - EANC * ( 1.0 - rT ) / Vtv ), 1.0 / NCN  );
	const auto IBEIP = IBEIP_TNOM * pow( pow(rT, XII) * exp ( - EAIC * ( 1.0 - rT ) / Vtv ), 1.0 / NCI  );
	const auto IBENP = IBENP_TNOM * pow( pow(rT, XIN) * exp ( - EANC * ( 1.0 - rT ) / Vtv ), 1.0 / NCN  );
	const auto IBCIP = IBCIP_TNOM * pow( pow(rT, XII) * exp ( - EAIS * ( 1.0 - rT ) / Vtv ), 1.0 / NCIP );
	const auto IBCNP = IBCNP_TNOM * pow( pow(rT, XIN) * exp ( - EANS * ( 1.0 - rT ) / Vtv ), 1.0 / NCNP );

//	Linear temperature mappings for NF/NR and AVC2
//		Note: this is an undesirable type of temperature mapping,
//		as it cannot be done "in-place" because it does not have the
//		properties P(T1->T2->T3)=P(T1->T3) and P(T1->T2->T1)=P(T1).
//		For this model it is best to always map from P_TNOM.

	const auto NF    =   NF_TNOM   * ( 1.0 + TNF   * ( Tdev - Tini ) );
	const auto NR    =   NR_TNOM   * ( 1.0 + TNF   * ( Tdev - Tini ) );
	const auto AVC2  = AVC2_TNOM   * ( 1.0 + TAVC  * ( Tdev - Tini ) );

//	Temperature mappings for built-in potentials

	const auto PE    = psibi ( PE_TNOM, EAIE, Vtv, rT );
	const auto PC    = psibi ( PC_TNOM, EAIC, Vtv, rT );
	const auto PS    = psibi ( PS_TNOM, EAIS, Vtv, rT );

//	zero-bias capacitance temperature mappings come directly from the
//	first-order theory for p-n junction capacitance as:

	const auto CJE   = CJE_TNOM  * pow( PE_TNOM / PE, ME);
	const auto CJC   = CJC_TNOM  * pow( PC_TNOM / PC, MC);
	const auto CJEP  = CJEP_TNOM * pow( PC_TNOM / PC, MC);
	const auto CJCP  = CJCP_TNOM * pow( PS_TNOM / PS, MS);

//	Temperature mappings for epi parameters

	const auto GAMM  = GAMM_TNOM * ( pow(rT, XIS) * exp ( - EA   * ( 1.0 - rT ) / Vtv ) );
	const auto VO    = VO_TNOM   * pow(rT, XVO);

//	End of temperature mappings.
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
		{	return 0.0;
		},
		[&](void)
		{	return Vrcx / RCX;
		}
	);
	const auto Irci_Kbci_Kbcx = if_(
		RCI <= 0.0,
		[](void)
		{	return std::make_tuple(0.0, 0.0, 0.0);
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
		{	return 0.0;
		},
		[&](void)
		{	return Vrbx / RBX;
		}
	);
	const auto Irbi = if_(
		RBI  <= 0.0,
		[](void)
		{	return 0.0;
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
		{	return 0.0;
		},
		[&](void)
		{	return Vre  / RE;
		}
	);
	const auto Irs = if_(
		RS   <= 0.0,
		[](void)
		{	return 0.0;
		},
		[&](void)
		{	return Vrs  / RS;
		}
	);
	const auto Irbp = if_(
		RBP <= 0.0,
		[](void)
		{	return 0.0;
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

		return std::make_tuple(
#define __create__(a, b, c) a
#define __COMMA__ ,
#include "outputs.h"
		);
//	End of equations.
	}
	static const char*const s_aNames[];
};
const char*const vbic::s_aNames[] = {
#define __create__(a, b, c) #a
#define __COMMA__ ,
#include "outputs.h"
};
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
}
int main(int, char**)
{	using namespace vbic95;
	vbic sI;
	std::array<double, std::size_t(enumNodes::NumberOfNodes)> sV({});
	std::cout << sI.calculate(sV) << "\n";
}
