#ifdef __DIODE__
	__create__(ID0, a0, NumberOfNodes) __COMMA__
	__create__(IR0, a0, b0) __COMMA__
	__create__(ID1, a1, NumberOfNodes) __COMMA__
	__create__(IR1, a1, b1)
#else
#ifdef SELF_HEATING
#ifdef EXCESS_PHASE
	__create__(Itxf, ci, ei) __COMMA__ //	Vrxf,delT	// forward   transport current
#else
	__create__(Itzf, ci, ei) __COMMA__ //		Vbei,Vbci,delT	// forward   transport current
#endif
	__create__(Itzr, ei, ci) __COMMA__ //		Vbei,Vbci,delT	// reverse   transport current
	__create__(Ibe, bi, ei) __COMMA__ //		Vbei,delT	// intrinsic b-e current
	__create__(Ibex, bx, ei) __COMMA__ //		Vbex,delT	// side      b-e current
	__create__(Ibc, bi, ci) __COMMA__ //		Vbci,delT	// intrinsic b-c current
	__create__(Igc, ci, bi) __COMMA__ //		Vbei,Vbci,delT	// c-b weak avalanche current
	__create__(Ircx, c, cx) __COMMA__ //		Vrcx,delT	// RCX  element
	__create__(Irci, cx, ci) __COMMA__ //		Vbci,Vrci,delT	// RCI  element
	__create__(Irbx, b, bx) __COMMA__ //		Vrbx,delT	// RBX  element
	__create__(Irbi, bx, bi) __COMMA__ //		Vrbi,Vbei,Vbci,delT	// RBI  element
	__create__(Ire, e, ei) __COMMA__ //		Vre,delT	// RE   element
	__create__(Irs, s, si) __COMMA__ //		Vrs,delT	// RS   element
	__create__(Iccp, bx, si) __COMMA__ //		Vbep,Vbcp,Vbci,delT	// parasitic transprt I
	__create__(Ibep, bx, bp) __COMMA__ //		Vbep,delT	// parasitic b-e current
	__create__(Ibcp, si, bp) __COMMA__ //		Vbcp,delT	// parasitic b-c current
	__create__(Irbp, cx, bp) __COMMA__ //		Vrbp,Vbep,Vbci,delT	// RBP element
	__create__(Ith, tl, dt) __COMMA__ //		Vbei,Vbci,Vrcx,Vrci,Vrbx,Vrbi,Vre,Vbep,Vbcp
				//Vrbp,Vrs,Vbex,delT
						// thermal power generation
	__create__(Irth, dt, tl) //		delT		// thermal resistance RTH
#else
#ifdef EXCESS_PHASE
	__create__(Itxf, ci, ei) __COMMA__ //		Vrxf		// forward   transport current
#else
	__create__(Itzf, ci, ei) __COMMA__ //		Vbei,Vbci	// forward   transport current
#endif
	__create__(Itzr, ei, ci) __COMMA__ //		Vbei,Vbci	// reverse   transport current
	__create__(Ibe, bi, ei) __COMMA__ //		Vbei		// intrinsic b-e current
	__create__(Ibex, bx, ei) __COMMA__ //		Vbex		// side      b-e current
	__create__(Ibc, bi, ci) __COMMA__ //		Vbci		// intrinsic b-c current
	__create__(Igc, ci, bi) __COMMA__ //		Vbei,Vbci	// c-b weak avalanche current
	__create__(Ircx, c, cx) __COMMA__ //		Vrcx		// RCX  element
	__create__(Irci, cx, ci) __COMMA__ //		Vbci,Vrci	// RCI  element
	__create__(Irbx, b, bx) __COMMA__ //		Vrbx		// RBX  element
	__create__(Irbi, bx, bi) __COMMA__ //	Vrbi,Vbei,Vbci	// RBI  element, it gets the
						//  Vbe/ci dependence through qb
	__create__(Ire, e, ei) __COMMA__ //	Vre		// RE   element
	__create__(Irs, s, si) __COMMA__ //	Vrs		// RS   element
	__create__(Iccp, bx, si) __COMMA__ //	Vbep,Vbcp,Vbci	// parasitic transport current
	__create__(Ibep, bx, bp) __COMMA__ //	Vbep		// parasitic b-e current
	__create__(Ibcp, si, bp) __COMMA__ //	Vbcp		// parasitic b-c current
	__create__(Irbp, cx, bp) //	Vrbp,Vbep,Vbci	// RBP element, it gets the
						//  Vbep and Vbci dependence
						//  through qbp
#endif
#ifdef EXCESS_PHASE
	__COMMA__
	__create__(Itzf, NumberOfNodes, xf1) __COMMA__ //, xf1 	Vbei,Vbci	// fwd transport zero phase
	__create__(Itxf, xf2,  NumberOfNodes) __COMMA__ // 	Vrxf		// fwd transport excess-phase
	__create__(Ilxf, xf1,  xf2)  // 			// excess-phase inductor, which
						// is the formulation variable
						// for this element and so is
						// not dependent on any branch
						// voltages
#endif
#endif
#undef __create__
#undef __COMMA__