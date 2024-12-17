# Two Classes for Automatic Differentiation

**Author**: Peter Foelsche
**Date**: October 2024
**Location**: Austin, TX, USA
**Email**: [peter_foelsche@outlook.com](mailto:peter_foelsche@outlook.com)

## Introduction

This document describes the implementation and usage of two classes designed for automatic differentiation leveraging dual numbers. These classes are designed for high performance, making use of sparse representations and template metaprogramming.

## Classes Overview

### **taylor::ctaylor**
- **Header File**: `ctaylor.h`
- **Description**: Implements a dual number truncated Taylor series class.
- **Template Parameter**: Max order of derivatives, should be larger than 1. For MAX=1, use `cjacobian.h`.
- **Mixing Instances**: Do not mix instances with different MAX parameters.
- **Functionality**: This class implements calculations with polynomial coefficients.

### **jacobian::cjacobian**
- **Header File**: `cjacobian.h`
- **Description**: Implements a dual number class for automatic derivation of max order 1. Much simpler than `ctaylor`.

## Implementation Details

Both implementations are sparse, carrying and calculating only potentially nonzero derivatives. Variables cannot be reused arbitrarily due to the type changing when different independent variables or derivative orders are involved. If-statements must be implemented using the `if_()` function, which supports returning tuples. The `auto` keyword should be used to allow the compiler to determine the correct type.

## Compile Time and Requirements

- **Compile Time**: Under Visual C++ (compared to g++ or clang++), especially for ctaylor and MAX > 1, compile time tends to be much Longer or even not finishing. (https://developercommunity.visualstudio.com/t/Compile-time-for-project-using-boost::mp/10760473). This does not apply for cjacobian or ctaylor with MAX=1.
- **C++ Standard**: Requires C++14.
- **Dependencies**: Requires `boost::mp11` (boost_1_86_0).
- **Note**: No double-cast operator to facilitate compiler errors for unimplemented functions.

## Test Cases

1. **ctaylor.cpp**
   - **Output**: `ctaylor.exe`
   - **Header**: `ctaylor.h`
   - **Description**: Reads arguments from the command line to prevent g++ from incorporating compile-time results into the executable. Used for checking correctness together with maxima.txt. Consult this example to learn how to create independent variables.
   - ** builds on Visual C++

2. **cjacobian.cpp**
   - **Output**: `cjacobian.exe`
   - **Header**: `cjacobian.h`
   - **Description**: Reads arguments from the command line. A very primitive implementation. Used for checking correctness together with maxima.txt. Consult this example to learn how to create independent variables.
   - ** builds on Visual C++

3. **VBIC95Jac.exe and VBIC95Taylor.exe**
   - **Headers**: `VBIC95/VBIC95.cpp` or `VBIC95Jac/VBIC95Jac.cpp`
   - **Description**: Implements a DC solver for a single transistor using VBIC95. Parameters are read from `VBIC95/PARS`.
   - ** builds on Visual C++
   - ** When attempting to use Halley's method in `VBIC95Taylor.exe` (increase `MAX` to 2 in `VBIC95/VBIC95.cpp`) Visual C++ hangs (filed a bug).

4. **BLACK_SCHOLES**
   - **Outputs**: `black_scholes.exe` and `black_scholes_orig.exe`
   - **Headers**: `ctaylor.h` or `boost/autodiff`
   - **Description**: Implements performance measurement (using an optional loop) and test for correctness (comparison with boost::autodiff).
   - ** Visual C++ Compiler hangs when trying to build (filed a bug).

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
        mp_list<mp_size_t<3>, mp_size_t<2>>
    > EXAMPLE; // 2nd derivative of independent variable with enum 3
    ```

  - For `cjacobian`:
    ```cpp
    template<std::size_t I>
    double cjacobian::getDer(const mp_size_t<I>&) const;
    ```

- **Handling Loops**: If `if_()` is insufficient for creating types in loops, refer to its implementation and use a similar approach manually. `taylor::common_type` and `jacobian::common_type` are used to merge `std::tuple` containing `taylor::ctaylor` or `jacobian::cjacobian`.

## Careful!

When using `ctaylor`, the maximum number of coefficients for a maximum order $N$ with $M$ independent variables is equal to $\binom{N+M}{M}$, which expands to $\frac{(N+M)!}{M!N!}$. Thus, one should be careful when increasing either the maximum order $N$ or the number of independent variables $M$, as the cost of storage, CPU-time, and compile-time increases quite dramatically.

I'm aware of the compile-time problems when using `ctaylor`—especially when using Visual C++—and I'm thinking about ways to make this less of a problem. However, I believe the resulting unrivaled performance is worth the trouble.

## Conclusion

These classes provide efficient and flexible implementations for automatic differentiation and truncated Taylor series calculations. By leveraging sparse representations and template metaprogramming, they ensure high performance and accuracy. Proper usage and understanding of the types and dependencies are crucial for integrating these classes into your projects effectively.
