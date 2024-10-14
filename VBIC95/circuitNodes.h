#ifdef __DIODE__
__create__(a0) __COMMA__
__create__(b0) __COMMA__
__create__(IVB0) __COMMA__
__create__(a1) __COMMA__
__create__(b1) __COMMA__
__create__(IVB1) __COMMA__
#else
__create__(c) __COMMA__
__create__(b) __COMMA__
__create__(e) __COMMA__
__create__(s) __COMMA__
#ifdef SELF_HEATING
__create__(dt) __COMMA__
//	__create__(tl) __COMMA__
#endif
__create__(cx) __COMMA__
__create__(ci) __COMMA__
__create__(bx) __COMMA__
__create__(bi) __COMMA__
__create__(ei) __COMMA__
__create__(si) __COMMA__
__create__(bp) __COMMA__
#ifdef EXCESS_PHASE
__create__(xf1) __COMMA__
__create__(xf2) __COMMA__
#endif
__create__(IVC) __COMMA__
__create__(IVB) __COMMA__
__create__(IVE) __COMMA__
__create__(IVS) __COMMA__
#endif
__create__(NumberOfNodes)
#undef __create__
#undef __COMMA__
