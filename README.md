# Two Classes for Automatic Differentiation

**Author**: Peter Foelsche |
**Date**: October 2024..January 2025 |
**Location**: Austin, TX, USA |
**Email**: [peter_foelsche@outlook.com](mailto:peter_foelsche@outlook.com)

## Introduction

This document describes the implementation and usage of two classes designed for automatic differentiation leveraging dual numbers. These classes are designed for high performance, making use of sparse representations and template metaprogramming.

## News

I attempted to get CI working, but abandoned the idea as it seems to be undocumented and unmanageable.
Anyway, this is a header-only library and what needs to be built here, are the various examples and the test executable.

Added CMakeLists.txt

Implemented the following functions:
- exp2()
- log2()
- isinf() -- returning if any value or derivative is infinite
- expm1()
- log1p()

Implemented tgamma() and polygamma().

Fixed undefined symbols which occured only for debug builds due to missing out-of-class definition of various static value objects.

I simplified ctaylor::chainRule(), reducing the number of steps.

For the `vbic*.exe`, I added support for asynchronous exception handling for floating-point exceptions, which is now finally available on Linux.

I added some automatic regression test (test/test.exe) using boost::test.

I changed the Makefile to enable the user to override `USER_CXXFLAGS=-march=native` with his own setting and I added a compiler flag `-fno-stack-protector` as this project is interested in performance but not in avoiding hacker attacks.

There is another compiler flag of interest not currently used (`-ffast-math`) which causes some testcases to run dramatically faster and others to be slightly slower.

I implemented and documented (in cjacobian.cpp and ctaylor.cpp) a way to perform chain-rule optimization.

jacobian::hypot() and potentially other nonlinear functions did not compile using g++ (why didn't anybody complain?!).
To fix this, I renamed all the nonlinear helper functions returning a std::pair to use a trailing underscore in the name, so that they wouldn't collide with the ordinary nonlinear function expecting and returning a cjacobian.

I fixed some compile time issue which prevented some projects using ctaylor.h to finish compiling on Visual C++.

## Classes Overview

### **taylor::ctaylor**
- **Header File**: `ctaylor.h`
- **Description**: Implements a dual number truncated Taylor series class for calculation of higher order derivatives.
- **Template Parameter**: Max order of derivatives, should be larger than 1. For MAX=1, use `cjacobian.h`.
- **Mixing Instances**: Do not attempt to mix instances with different MAX parameters.
- **Functionality**: This class implements calculations with polynomial coefficients.

### **jacobian::cjacobian**
- **Header File**: `cjacobian.h`
- **Description**: Implements a dual number class for automatic derivation of max order 1. Much simpler than `ctaylor`.

## bulding the regression test and the examples

This is not strictly necessary, as using these classes can be done simply by including them.
There are 3 ways to build the examples and the regression test:

- **Visual C++ 2022**: load `ctaylor.sln` after having set `%BOOST_ROOT%` and potentially having created a symbolic link called `%BOOST_ROOOT%\include` pointing to `%BOOST_ROOT%` -- see LINUX-make
- **LINUX-make**: after having set the environment variable `BOOST_ROOT` to the location of the boost directory. I added a subdirectory `$(BOOST_ROOT)/include`, which in case of you have not let the build tool of boost (`b2.exe`) install the boost build, can simply be a symbolic link to `$(BOOST_ROOT)` itself -- even on Windows.
```
rem on windows as administrator
cd %BOOST_ROOT%
mklink /D include .
```
```
#on LINUX
cd $(BOOST_ROOT)
ln -s . include
```
build:
```
cd ctaylor
make -j$(nproc)
cd test
./test.exe
```
- **cmake**:
```
cd ctaylor
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --verbose -j$(nproc)
ctest
```
## Implementation Details

Both implementations are sparse, carrying and calculating only potentially nonzero derivatives. Variables cannot be reused arbitrarily due to the type changing when different independent variables or derivative orders are involved. If-statements must be implemented using the `if_()` function, which supports returning tuples. The `auto` keyword should be used to allow the compiler to determine the correct type. Don't specify the the full `ctaylor`/`cjacobian` type by hand, but use `decltype()` or `jacobian::common_type` or `taylor::common_type` or best use `auto`.

## Compile Time and Requirements

- **Compile Time**: Under Visual C++ (compared to g++ or clang++), especially for ctaylor and MAX > 1, compile time tends to be much longer or even infinite. I filed a <a href="https://developercommunity.visualstudio.com/t/Compile-time-for-project-using-boost::mp/10760473">bug</a> regarding this problem. This does not apply for cjacobian or ctaylor with MAX=1 or any of the testcases provided here.
- **C++ Standard**: Requires C++14.
- **Dependencies**: Requires `boost::mp11` (boost_1_86_0). Use the environment variable BOOST_ROOT to specify the location of the boost include files.
- **Note**: No double-cast operator to facilitate compiler errors for unimplemented functions.

## Test Cases

1. **ctaylor.cpp**
   - **Output**: `ctaylor.exe`
   - **Header**: `ctaylor.h`
   - **Description**: Reads arguments from the command line to prevent g++ from incorporating compile-time results into the executable. Used for checking correctness together with maxima.txt which represents an input file for <a href="https://maxima.sourceforge.io/">maxima</a>. Consult this example to learn how to create independent variables. Consult this example to learn how to perform chain-rule optimization.
   - **Visual C++**: ok

2. **cjacobian.cpp**
   - **Output**: `cjacobian.exe`
   - **Header**: `cjacobian.h`
   - **Description**: Reads arguments from the command line. A very primitive implementation. Used for checking correctness together with maxima.txt. Consult this example to learn how to create independent variables. Consult this example to learn how to perform chain-rule optimization.
   - **Visual C++**: ok

3. **VBIC95/VBIC95.cpp** and **VBIC95Jac/VBIC95Jac.cpp**
   - **Outputs**: `vbic95Jac.exe` and `vbic95Taylor.exe`
   - **Headers**: `cjacobian.h` or `ctaylor.h`
   - **Description**: Implements a DC solver for a single transistor using VBIC95. Parameters are read from `VBIC95/PARS`.
   - **Visual C++**: ok. When attempting to use Halley's method in `vbic95Taylor.exe` (increase `MAX` to 2 in `VBIC95/VBIC95.cpp`) Visual C++ takes a some time (4min:36s on my 12year old Dell M4700).
   - **g++-13**: 0min:13s (vbic95Taylor.exe for MAX=2)

4. **BLACK_SCHOLES**
   - **Outputs**: `black_scholes.exe` and `black_scholes_orig.exe`
   - **Headers**: `ctaylor.h` or `boost/autodiff`
   - **Description**: Implements performance measurement (using an optional loop) and test for correctness (comparison with boost::autodiff).
   - **Visual C++**: `black_scholes.exe` takes 9min:24s on my 12 year old Dell M4700.
   - **g++-13**: 0min:14s
   - **performance test ctaylor**: `time black_scholes.exe 0.02` shows 18s
   - **performance test autodiff**: `time ./black_scholes_orig.exe 0.02` shows 260s. This slow performance can partly be explained by the fact that `boost::autodiff` calculates many more derivatives. It does not consider the sum of all derivative orders for a term. Instead, the order with respect to any independent variable can be set independently of the others. As a result, the maximum order using `boost::autodiff` is the sum of all maximum orders of each independent variable. The method used by `boost::autodiff` for setting orders separately for every independent variable makes it ineffective for any use case with more than a single independent variable and a maximum order larger than one.

5. **logistic_regression**
   - **Outputs**: `logistic_regression.exe`
   - **Headers**: `cjacobian.h`
   - **Description**: Supposed to prove that forward-mode AD is not necessarily slower than reverse-mode AD but quite the opposite when applying the chain rule properly. Compare to [reverse-mode AD](https://github.com/ExcessPhase/reverse_mode_automatic_differentiation).
   - **Visual C++**: no problem
   - **g++-13**: no problem

6. **boost test regression test cases**
   - **Outputs**: `test/test.exe`
   - **Headers**: `cjacobian.h` and `ctaylor.h`
   - **Description**: .Compares results with `test/data*.txt` files which have been precalculated using maxima and `test/maxima*.txt`. Should be executed with the current directory equal to `test` in order to find `data*.txt`!
   - **Visual C++**: no problem
   - **g++-13**: no problem

## Accessing Results

- **0th Derivative (Value)**: Use `value(source)` in both `ctaylor.h` and `cjacobian.h`.
- **Specific Derivatives**:
  - For `ctaylor`:
    ```cpp
    template<typename LIST>
    double ctaylor::getDer(const LIST&) const;
    ```
    Example:
    ```cpp
    typedef mp_list<
        pair<mp_size_t<3>, mp_size_t<2>>
    > EXAMPLE; // 2nd derivative of independent variable with enum 3
    ```

  - For `cjacobian`:
    ```cpp
    template<std::size_t I>
    double cjacobian::getDer(const mp_size_t<I>&) const;
    ```

## Chain-Rule optimization

Imagine `f(g(x, y, z, t))`, with `x`, `y`, `z`, and `t` being independent variables. Let's also imagine that `f(g)` is non-trivial and would benefit from calculating only one derivative with respect to `g()`, instead of four derivatives with respect to `x`, `y`, `z`, and `t`.

In this case, one could store the result of `g(x, y, z, t)` as a variable and create a new independent variable `gx`. Then, evaluate `f(gx)` and for this value restore the correct derivatives with respect to `x`, `y`, `z`, and `t`.

When doing so, one should use a unique enumeration value that does not overlap with any other enumeration values in use. One can use the new independent variable instead of the old ones, and in the final result, the chain rule should be applied to go back to the original variables.
```
	/// original wasteful code
const auto fx = f(g(x, y, z, t));
```
```
	/// code leveraging chain-rule optimization
const auto gx = g(x, y, z, t);
const auto gx1 = gx.convert2Independent(mp_size_t<ENUM_G>());
const auto fx = f(gx1).chainRule(gx, mp_size_t<ENUM_G>());
```

## Handling Loops

If `if_()` is insufficient for creating types in loops, refer to its implementation and use a similar approach manually. `taylor::common_type` and `jacobian::common_type` are used to merge `std::tuple` containing `taylor::ctaylor` or `jacobian::cjacobian`.

## Careful!

When using `ctaylor`, the maximum number of coefficients for a maximum order $N$ with $M$ independent variables is equal to $\binom{N+M}{M}$, which expands to $\frac{(N+M)!}{M!N!}$. Thus, one should be careful when increasing either the maximum order $N$ or the number of independent variables $M$, as the cost of storage, CPU-time, and compile-time increases quite dramatically.

I'm aware of the compile-time problems when using `ctaylor`—especially when using Visual C++—and I'm thinking about ways to make this less of a problem. However, I believe the resulting unrivaled performance is worth the trouble.

## Conclusion

These classes provide efficient and flexible implementations for automatic differentiation and truncated Taylor series calculations. By leveraging sparse representations and template metaprogramming, they ensure high performance and accuracy. Proper usage and understanding of the types and dependencies are crucial for integrating these classes into your projects effectively.
