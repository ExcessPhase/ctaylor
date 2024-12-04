# Dual Number Truncated Taylor Series and Automatic Differentiation Classes

**Author**: Peter Foelsche
**Date**: October 2024
**Location**: Austin, TX, USA
**Email**: [peter_foelsche@outlook.com](mailto:peter_foelsche@outlook.com)

## Introduction

This document describes the implementation and usage of two classes designed for automatic differentiation and dual number truncated Taylor series calculations. These classes are designed for high performance and accuracy, making use of sparse representations and template metaprogramming.

## Classes Overview

### **taylor::ctaylor**
- **Header File**: `ctaylor.h`
- **Description**: Implements a dual number truncated Taylor series class.
- **Template Parameter**: Max order of derivatives, must be larger than 1. For MAX=1, use `cjacobian.h`.
- **Mixing Instances**: Do not mix instances with different MAX parameters.
- **Functionality**: This class implements calculations with polynomial coefficients.

### **jacobian::cjacobian**
- **Header File**: `cjacobian.h`
- **Description**: Implements a dual number class for automatic derivation of max order 1. Much simpler than `ctaylor`.

## Implementation Details

Both implementations are sparse, carrying and calculating only potentially nonzero derivatives. Variables cannot be reused arbitrarily due to type changes with independent variables or derivative orders. If-statements must be implemented using the `if_()` function, which supports tuples. The `auto` keyword should be used to allow the compiler to determine the correct type.

## Compile Time and Requirements

- **Compile Time**: Longer under Visual C++ 2022 compared to g++.
- **C++ Standard**: Requires C++14.
- **Dependencies**: Requires `boost::mp11` (boost_1_86_0).
- **Note**: No double-cast operator to facilitate compiler errors for unimplemented functions.

## Test Cases

1. **ctaylor.cpp**
   - **Output**: `ctaylor.exe`
   - **Header**: `ctaylor.h`
   - **Description**: Reads arguments from the command line to prevent g++ from incorporating compile-time results into the executable.

2. **cjacobian.cpp**
   - **Output**: `cjacobian.exe`
   - **Header**: `cjacobian.h`
   - **Description**: Reads arguments from the command line. A very primitive implementation.

3. **VBIC95Jac.exe and VBIC95Taylor.exe**
   - **Headers**: `VBIC95/VBIC95.cpp` or `VBIC95Jac/VBIC95Jac.cpp`
   - **Description**: Implements a DC solver for a single transistor using VBIC95. Parameters are read from `VBIC95/PARS`.

4. **BLACK_SCHOLES**
   - **Outputs**: `black_scholes.exe` and `black_scholes_orig.exe`
   - **Headers**: `ctaylor.h` or `boost/autodiff`
   - **Description**: Implements performance measurement using an optional loop.

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

## Conclusion

These classes provide efficient and flexible implementations for automatic differentiation and truncated Taylor series calculations. By leveraging sparse representations and template metaprogramming, they ensure high performance and accuracy. Proper usage and understanding of the types and dependencies are crucial for integrating these classes into your projects effectively.
