published under MIT license

Author: Peter Foelsche

October-1th 2024
Austin, TX, USA

email:	peter_foelsche@outlook.com

A sparse, dual number implementation for calculating not just the 1th order of derivatives.
Using this class to calculate derivative of max. 1th order, would constitute a waste,
as a dual number implementation for 1th order derivative is much simpler and cheaper.
Refer to ctaylor.cpp for a example usage.
Compile time under Visual C++ 2022 tends to be much longer than using g++.
Requires C++14.
Requires boost::mp11 (boost_1_86_0).

There are two examples:

1) 
ctaylor.cpp
Reads arguments from the commandline as otherwise the g++ compiler simply incorporates the at compile-time calculated results into the executable.
2) 
VBIC95.cpp -- simulates a single transistor using Halley's method (used to use Newton's method).
Reads parameters from file named PARS.
Still not fully tested.

The type of the dual numbers object depends on how many independent variables are involved and their cross-products and order.
So it is a good idea to avoid declaring variables in a different way than letting the compiler decide the type by using auto.
Reusing a variable for different purposes, is likely to fail (see the original bsim3 code).
For joining two different types, I created the taylor::if_() function, which operates like an terniary operator -- with double and ctaylor and tuples of such.
I think I've yet to implement an assignment operator -- which would look similar like the copy-constructor.

10/10/2024:
implemented assignment operators and cbrt()
The factorization algorithm in LUFAC does not attempt to achieve performance or perfect ordering.
In fact, it is the reason why VBIC fails at higher currents.
VBIC simulator can be changed to include self-heating by defining SELF_HEATING and excess-phase by defining EXCESS_PHASE
Copy constructors and assignment operators between different types only compile, if the target type contains all the values of the source type.
Means that it is impossible to forget by mistake derivatives.

In order to access the value (0th derivative) use 
friend double value(const ctaylor&).
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
