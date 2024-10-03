//	This section defines mappings of temperature dependent parameters.
//	Note that with SELF-HEATING these mappings must be done at each
//	bias calculation, and cannot be relegated to preprocessing, because
//	the "voltage" for the local temperature rise delT must be taken
//	into account.

#ifdef SELF_HEATING
	__CONST_AUTO__ Tdev  = TAMB + TABS + delT;	// device temperature (K), with delT
#else
	__CONST_AUTO__ Tdev  = TAMB + TABS;		// device temperature (K)
#endif
	__CONST_AUTO__ Tini  = TNOM + TABS;		// TNOM in K

	__CONST_AUTO__ Vtv   = KB * Tdev / QQ;		// thermal voltage
	__CONST_AUTO__ rT    = Tdev / Tini;		// ratio Tdev/TNOM in K

	__CONST_AUTO__ RCX   = RCX_TNOM * pow(rT, XRC);	// temperature mapping
	__CONST_AUTO__ RCI   = RCI_TNOM * pow(rT, XRC);	// of resistances is
	__CONST_AUTO__ RBX   = RBX_TNOM * pow(rT, XRB);	// according to mobility
	__CONST_AUTO__ RBI   = RBI_TNOM * pow(rT, XRB);	// variation, which
	__CONST_AUTO__ RE    = RE_TNOM  * pow(rT, XRE);	// should be same in
	__CONST_AUTO__ RS    = RS_TNOM  * pow(rT, XRS);	// emit/coll & base/subs
	__CONST_AUTO__ RBP   = RBP_TNOM * pow(rT, XRC);	// Note: RBP is in coll

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

	__CONST_AUTO__ IS    = IS_TNOM    * pow( pow(rT, XIS) * exp ( - EA   * ( 1.0 - rT ) / Vtv ), 1.0 / NF_TNOM );
	__CONST_AUTO__ ISP   = ISP_TNOM   * pow( pow(rT, XIS) * exp ( - EA   * ( 1.0 - rT ) / Vtv ), 1.0 / NFP  );
	__CONST_AUTO__ IBEI  = IBEI_TNOM  * pow( pow(rT, XII) * exp ( - EAIE * ( 1.0 - rT ) / Vtv ), 1.0 / NEI  );
	__CONST_AUTO__ IBEN  = IBEN_TNOM  * pow( pow(rT, XIN) * exp ( - EANE * ( 1.0 - rT ) / Vtv ), 1.0 / NEN  );
	__CONST_AUTO__ IBCI  = IBCI_TNOM  * pow( pow(rT, XII) * exp ( - EAIC * ( 1.0 - rT ) / Vtv ), 1.0 / NCI  );
	__CONST_AUTO__ IBCN  = IBCN_TNOM  * pow( pow(rT, XIN) * exp ( - EANC * ( 1.0 - rT ) / Vtv ), 1.0 / NCN  );
	__CONST_AUTO__ IBEIP = IBEIP_TNOM * pow( pow(rT, XII) * exp ( - EAIC * ( 1.0 - rT ) / Vtv ), 1.0 / NCI  );
	__CONST_AUTO__ IBENP = IBENP_TNOM * pow( pow(rT, XIN) * exp ( - EANC * ( 1.0 - rT ) / Vtv ), 1.0 / NCN  );
	__CONST_AUTO__ IBCIP = IBCIP_TNOM * pow( pow(rT, XII) * exp ( - EAIS * ( 1.0 - rT ) / Vtv ), 1.0 / NCIP );
	__CONST_AUTO__ IBCNP = IBCNP_TNOM * pow( pow(rT, XIN) * exp ( - EANS * ( 1.0 - rT ) / Vtv ), 1.0 / NCNP );

//	Linear temperature mappings for NF/NR and AVC2
//		Note: this is an undesirable type of temperature mapping,
//		as it cannot be done "in-place" because it does not have the
//		properties P(T1->T2->T3)=P(T1->T3) and P(T1->T2->T1)=P(T1).
//		For this model it is best to always map from P_TNOM.

	__CONST_AUTO__ NF    =   NF_TNOM   * ( 1.0 + TNF   * ( Tdev - Tini ) );
	__CONST_AUTO__ NR    =   NR_TNOM   * ( 1.0 + TNF   * ( Tdev - Tini ) );
	__CONST_AUTO__ AVC2  = AVC2_TNOM   * ( 1.0 + TAVC  * ( Tdev - Tini ) );

//	Temperature mappings for built-in potentials

	__CONST_AUTO__ PE    = psibi ( PE_TNOM, EAIE, Vtv, rT );
	__CONST_AUTO__ PC    = psibi ( PC_TNOM, EAIC, Vtv, rT );
	__CONST_AUTO__ PS    = psibi ( PS_TNOM, EAIS, Vtv, rT );

//	zero-bias capacitance temperature mappings come directly from the
//	first-order theory for p-n junction capacitance as:

	__CONST_AUTO__ CJE   = CJE_TNOM  * pow( PE_TNOM / PE, ME);
	__CONST_AUTO__ CJC   = CJC_TNOM  * pow( PC_TNOM / PC, MC);
	__CONST_AUTO__ CJEP  = CJEP_TNOM * pow( PC_TNOM / PC, MC);
	__CONST_AUTO__ CJCP  = CJCP_TNOM * pow( PS_TNOM / PS, MS);

//	Temperature mappings for epi parameters

	__CONST_AUTO__ GAMM  = GAMM_TNOM * ( pow(rT, XIS) * exp ( - EA   * ( 1.0 - rT ) / Vtv ) );
	__CONST_AUTO__ VO    = VO_TNOM   * pow(rT, XVO);

//	End of temperature mappings.
#undef __CONST_AUTO__ 
