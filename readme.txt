published under MIT license

Author: Peter Foelsche

October-1th 2024
Austin, TX, USA

email:	peter_foelsche@outlook.com

There are two classes:
	1) taylor/ctaylor in ctaylor.h
		Implements a dual number truncated taylor series class.
		Max order of derivatives is a template parameter and should be larger than 1.
		For MAX=1 use cjacobian.h
	2) jacobian/cjacobian in cjacobian.h
		Implements a dual number class for automatic derivation of max. order 1.

Both implementations are sparse. They carry & calculate only these derivatives, which are actually nonzero.
This means, that variables cannot be reused at will because the types changes with either the involved independent variables (cjacobian.h)
or even the order and cross derivatives involved (ctaylor.h).
Also if there are multiple branches (if-statements).
They have to be realized using the if_() function, which also works for tuples.
See the vbic85 examples.
That's why it is a good idea to use the auto keyword to let the compiler select the correct type.

Compile time under Visual C++ 2022 tends to be much longer than using g++.
Requires C++14.
Requires boost::mp11 (boost_1_86_0).
There is intentionally no double-cast operator in order to facilitate a compiler error if an unimplemented function is being used.

There are multiple testcases:

1) ctaylor.cpp ctaylor.exe uses ctaylor.h
Reads arguments from the commandline as otherwise the g++ compiler simply incorporates the at compile-time calculated results into the executable.
2)
cjacobian.cpp yields cjacobian.exe using cjacobian.h
Reads arguments from the commandline.
Very primitive.
3)
vbic95Jac.exe
vbic95Taylor.exe
uses
	VBIC95/VBIC95.cpp
or
	VBIC95Jac/VBIC95Jac.cpp

Implements a DC solver for a single transistor using VBIC95.
Reads parameters from VBIC95/PARS which needs to be passed on the commandline.
Either using cjacobian.h or ctaylor.h with MAX=1.
vbic95Taylor.exe implements Halley's method which defaults to Newton's method in case of MAX=1 which causes the Hessian to be zero.
It does not converge for certain bias points when using MAX=2 (Halley's method). Potentially the scaling of delta-x does screw up Halley's method.
But even without scaling of delta-X it does not converge for Halley's method.

4)
BLACK_SCHOLES
yields black_scholes.exe
using ctaylor.h
Example code from the boost library. If it would be started in a loop without printouts, it would show a dramatic performance improvement compared to boost::autodiff.
This example is extracting values from complicated expressions involving ctaylor variables.
This constitutes unnecessary calculation of derivatives just to increase entropy of the universe.
I don't know if the compiler is smart enough to avoid this.
I was unable to compile this testcase using Visual C++ 2022.

Access to results:
For accessing the 0th derivative (value) use the value(source) function which is implemented in both, ctaylor.h and cjacobian.h.
In order to access a particular derivative, use

template<typename LIST>
double ctaylor::getDer(const LIST&) const;

LIST being an argument containing the derivative to access, e.g.

typedef mp_list<
	mp_list<
		mp_size_t<3>,
		mp_size_t<2>
	>
> EXAMPLE;

for 2nd derivative of independent variable with enum 3.

typedef mp_list<
	mp_list<
		mp_size_t<0>,
		mp_size_t<1>
	>,
	mp_list<
		mp_size_t<1>,
		mp_size_t<1>
	>
> EXAMPLE;

for d2/dx0/dx1!

template<std::size_t I>
double cjacobian::getDer(const mp_size_t<I>&) const

Pass the integer identifying the independent variable as a type instance.
