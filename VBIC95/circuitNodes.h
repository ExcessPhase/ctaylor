__create2__(c, -2) __COMMA__
__create2__(b, -1) __COMMA__
//__create__(e) __COMMA__
//__create__(s) __COMMA__
__create__(cx) __COMMA__
__create__(ci) __COMMA__
__create__(bx) __COMMA__
__create__(bi) __COMMA__
__create__(ei) __COMMA__
__create__(si) __COMMA__
__create__(bp) __COMMA__
#ifdef SELF_HEATING
__create__(dt) __COMMA__
//	__create__(tl) __COMMA__
#endif
#ifdef EXCESS_PHASE
__create__(xf1) __COMMA__
__create__(xf2) __COMMA__
#endif
__create__(NumberOfNodes)
#undef __create__
#undef __create2__
#undef __COMMA__
