#ifdef EXCESS_PHASE
__create__(Itxf) __COMMA__
#else
__create__(Itzf) __COMMA__
#endif
#ifdef SELF_HEATING
__create__(Ith) __COMMA__
__create__(Irth) __COMMA__
#endif
__create__(Itzr) __COMMA__
__create__(Ibe) __COMMA__
__create__(Ibex) __COMMA__
__create__(Ibc) __COMMA__
__create__(Igc) __COMMA__
__create__(Ircx) __COMMA__
__create__(Irci) __COMMA__
__create__(Irbx) __COMMA__
__create__(Irbi) __COMMA__
__create__(Ire) __COMMA__
__create__(Irs) __COMMA__
__create__(Iccp) __COMMA__
__create__(Ibep) __COMMA__
__create__(Ibcp) __COMMA__
__create__(Irbp)
#undef __create__
#undef __COMMA__
