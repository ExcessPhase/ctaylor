//	This section defines mappings of temperature dependent parameters.
//	Note that with SELF-HEATING these mappings must be done at each
//	bias calculation, and cannot be relegated to preprocessing, because
//	the "voltage" for the local temperature rise delT must be taken
//	into account.

#ifdef SELF_HEATING
	__create2__(Tdev, TAMB + TABS + delT) __COMMA__	// device temperature (K), with delT
#else
	__create2__(Tdev, TAMB + TABS) __COMMA__		// device temperature (K)
#endif
	__create2__(Tini, TNOM + TABS) __COMMA__		// TNOM in K

	__create2__(Vtv, KB * Tdev / QQ) __COMMA__		// thermal voltage
	__create2__(rT, Tdev / Tini) __COMMA__		// ratio Tdev/TNOM in K

	__create__(RCX, RCX_TNOM * pow(rT, XRC)) __COMMA__	// temperature mapping
	__create__(RCI, RCI_TNOM * pow(rT, XRC)) __COMMA__	// of resistances is
	__create__(RBX, RBX_TNOM * pow(rT, XRB)) __COMMA__	// according to mobility
	__create__(RBI, RBI_TNOM * pow(rT, XRB)) __COMMA__	// variation, which
	__create__(RE, RE_TNOM  * pow(rT, XRE)) __COMMA__	// should be same in
	__create__(RS, RS_TNOM  * pow(rT, XRS)) __COMMA__	// emit/coll & base/subs
	__create__(RBP, RBP_TNOM * pow(rT, XRC)) __COMMA__	// Note: RBP is in coll

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

	__create__(IS, IS_TNOM    * pow( pow(rT, XIS) * exp ( - EA   * ( 1.0 - rT ) / Vtv ), 1.0 / NF_TNOM )) __COMMA__
	__create__(ISP, ISP_TNOM   * pow( pow(rT, XIS) * exp ( - EA   * ( 1.0 - rT ) / Vtv ), 1.0 / NFP  )) __COMMA__
	__create__(IBEI, IBEI_TNOM  * pow( pow(rT, XII) * exp ( - EAIE * ( 1.0 - rT ) / Vtv ), 1.0 / NEI  )) __COMMA__
	__create__(IBEN, IBEN_TNOM  * pow( pow(rT, XIN) * exp ( - EANE * ( 1.0 - rT ) / Vtv ), 1.0 / NEN  )) __COMMA__
	__create__(IBCI, IBCI_TNOM  * pow( pow(rT, XII) * exp ( - EAIC * ( 1.0 - rT ) / Vtv ), 1.0 / NCI  )) __COMMA__
	__create__(IBCN, IBCN_TNOM  * pow( pow(rT, XIN) * exp ( - EANC * ( 1.0 - rT ) / Vtv ), 1.0 / NCN  )) __COMMA__
	__create__(IBEIP, IBEIP_TNOM * pow( pow(rT, XII) * exp ( - EAIC * ( 1.0 - rT ) / Vtv ), 1.0 / NCI  )) __COMMA__
	__create__(IBENP, IBENP_TNOM * pow( pow(rT, XIN) * exp ( - EANC * ( 1.0 - rT ) / Vtv ), 1.0 / NCN  )) __COMMA__
	__create__(IBCIP, IBCIP_TNOM * pow( pow(rT, XII) * exp ( - EAIS * ( 1.0 - rT ) / Vtv ), 1.0 / NCIP )) __COMMA__
	__create__(IBCNP, IBCNP_TNOM * pow( pow(rT, XIN) * exp ( - EANS * ( 1.0 - rT ) / Vtv ), 1.0 / NCNP )) __COMMA__

//	Linear temperature mappings for NF/NR and AVC2
//		Note: this is an undesirable type of temperature mapping,
//		as it cannot be done "in-place" because it does not have the
//		properties P(T1->T2->T3)=P(T1->T3) and P(T1->T2->T1)=P(T1).
//		For this model it is best to always map from P_TNOM.

	__create__(NF, NF_TNOM   * ( 1.0 + TNF   * ( Tdev - Tini ) )) __COMMA__
	__create__(NR, NR_TNOM   * ( 1.0 + TNF   * ( Tdev - Tini ) )) __COMMA__
	__create__(AVC2, AVC2_TNOM   * ( 1.0 + TAVC  * ( Tdev - Tini ) )) __COMMA__

//	Temperature mappings for built-in potentials

	__create__(PE, psibi ( PE_TNOM, EAIE, Vtv, rT )) __COMMA__
	__create__(PC, psibi ( PC_TNOM, EAIC, Vtv, rT )) __COMMA__
	__create__(PS, psibi ( PS_TNOM, EAIS, Vtv, rT )) __COMMA__

//	zero-bias capacitance temperature mappings come directly from the
//	first-order theory for p-n junction capacitance as:

	__create__(CJE, CJE_TNOM  * pow( PE_TNOM / PE, ME)) __COMMA__
	__create__(CJC, CJC_TNOM  * pow( PC_TNOM / PC, MC)) __COMMA__
	__create__(CJEP, CJEP_TNOM * pow( PC_TNOM / PC, MC)) __COMMA__
	__create__(CJCP, CJCP_TNOM * pow( PS_TNOM / PS, MS)) __COMMA__

//	Temperature mappings for epi parameters

	__create__(GAMM, GAMM_TNOM * ( pow(rT, XIS) * exp ( - EA   * ( 1.0 - rT ) / Vtv ) )) __COMMA__
	__create__(VO, VO_TNOM   * pow(rT, XVO)) __COMMA__

//	End of temperature mappings.
#undef __create__
#undef __create2__
#undef __COMMA__ 
