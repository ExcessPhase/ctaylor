published under MIT license

Author: Peter Foelsche

October-1th 2024
Austin, TX, USA

email:	peter_foelsche@outlook.com

A sparse, dual number implementation for calculating not just the 1th order of derivatives
Refer to ctaylor.cpp for a example usage
Compile time under Visual C++ 2022 tends to be much longer than using g++
Requires C++14
Requires boost::mp11 (boost_1_86_0)

There are two examples:

1) 
ctaylor.cpp
Reads arguments from the commandline as otherwise the g++ compiler simply incorporates the compile-time calculated results into the executable.
2) 
VBIC95.cpp -- simulates a single transistor using Halley's method. 
Reads parameters from file named PARS.
Still not fully tested.
