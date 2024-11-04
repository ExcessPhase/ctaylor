published under MIT license

Author: Peter Foelsche

October 2024
Austin, TX, USA

email:	peter_foelsche@outlook.com

There are two classes:
	1) taylor::ctaylor in ctaylor.h
		Implements a dual number truncated taylor series class.
		Max order of derivatives is a template parameter and should be larger than 1.
		For MAX=1 use cjacobian.h
		Don't attempt to mix instances of ctaylor with different MAX parameter.
		This class implements basically a calculation with polynomial coefficients.
	2) jacobian::cjacobian in cjacobian.h
		Implements a dual number class for automatic derivation of max. order 1.
		Much simpler than ctaylor

Both implementations are sparse. They carry & calculate only these derivatives, which could potentially be nonzero.
This means, that variables cannot be reused at will because the types changes with either the involved independent variables (cjacobian.h)
or even the order and cross derivatives involved (ctaylor.h).
Also if there are multiple branches (if-statements) the results of both arms might differ in type.
If-statements have to be realized using the if_() function, which also works for tuples.
See the vbic95 examples.
That's why it is a good idea to use the auto keyword to let the compiler select the correct type.

Compile time under Visual C++ 2022 tends to be much longer than using g++.
Requires C++14.
Requires boost::mp11 (boost_1_86_0).
There is intentionally no double-cast operator in order to facilitate a compiler error if an unimplemented function is being used.

There are multiple testcases:

1) ctaylor.cpp yields ctaylor.exe uses ctaylor.h
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
But even without scaling of delta-X it does not converge for Halley's method (for certain bias points).
There are as many versions of Halley's method in the internet as there are publications about it!

4)
BLACK_SCHOLES
yields black_scholes.exe
using ctaylor.h
Example code from the boost library. If it would be started in a loop without printouts, it would show a dramatic performance improvement compared to boost::autodiff.
Some expressions in this example are using ctaylor variables in non-trivial expressions and at the end only the value is being used by enclosing the entire expression in a call to value().
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

There might be one area in which if_() is not sufficient to create types, able to carry all assigned derivatives: loops

In this case one might have a look at the implementation of if_() and use this method by hand.
I implemented not std::common_type<> but taylor::common_type<>/jacobian::common_type as I wanted to merge std::tuple<> as well. Of course only std::tuple<> which contain a taylor::ctaylor<>/jacobian::cjacobian<>. Potentially I'm going to change this to use std::common_type<>.
