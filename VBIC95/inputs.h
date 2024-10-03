#ifdef SELF_HEATING
	__create__(delT, dt, tl) __COMMA__// voltage across RTH, local temperature rise
					// measured with respect to ground (ambient)
#endif
#ifdef EXCESS_PHASE
	__create__(Vcxf, xf1, NumberOfNodes) __COMMA__	// voltage across excess-phase capacitor
	__create__(Vrxf, xf2, NumberOfNodes) __COMMA__	// voltage across excess-phase resistor Itxf
#endif
	__create__(Vbe, b, e) __COMMA__
	__create__(Vbc, b, c) __COMMA__
	__create__(Vbei, bi, ei) __COMMA__	// intrinsic b-e voltage
	__create__(Vbex, bx, ei) __COMMA__	// side      b-e voltage
	__create__(Vbci, bi, ci) __COMMA__	// intrinsic b-c voltage
	__create__(Vrcx, c, cx) __COMMA__	// voltage across RCX
	__create__(Vrci, cx, ci) __COMMA__	// voltage across RCI
	__create__(Vrbx, b, bx) __COMMA__	// voltage across RBX
	__create__(Vrbi, bx, bi) __COMMA__	// voltage across RBI
	__create__(Vre, e, ei) __COMMA__	// voltage across RE
	__create__(Vrs, s, si) __COMMA__	// voltage across RS
	__create__(Vbep, bx, bp) __COMMA__	// parasitic b-e voltage (pnp polarity)
	__create__(Vbcp, si, bp) __COMMA__	// parasitic b-c voltage (pnp polarity)
	__create__(Vrbp, cx, bp) // voltage across RBP
#undef __create__
#undef __COMMA__