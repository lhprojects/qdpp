

# What is this library?

This library `qdpp` supports double-`double` (104bits) and quad-`double` (206bits) floating points.

This is library  is based on the QD library.

# Why is this library?

`qdpp` is based on the QD library. 

## So why `QD`?

`QD` supports double-`double` (104bits) and quad-`double` (206bits) floating points. The speed of double-`double` and quad-`double` are much faster than the normal multiple-precision library. Because

* Dynamic memory allocation is not need. the memory is allocated on stack.
* Easy to implement compiling time evaluation math functions.

## `qdpp` compared to `QD`

* Make the library header only for c++.
* Add `constexpr`  to math functions for `dd_real` and `qd_real`.\
* Improvements for some math and IO functions.
* `constexpr` math functions have been implemented for double.
* This is library is not 100% compatible to the QD library.

## Compute constants in compiling time
```c++
constinit dd_real dd_sqrt2 = sqrt(dd_real(2.));
constexpr dd_real dd_sqrt2_constexpr = sqrt(dd_real(2.));;

constinit double double_sqrt2 = fb::sqrt(2.);
constexpr double double_sqrt2_constexpr = fb::sqrt(2.);
constinit double double_sqrt2_sqr = double_sqrt2_constexpr*double_sqrt2_constexpr;

constinit double dd_sqrt2_sqr_to_double = to_double(dd_sqrt2_constexpr*dd_sqrt2_constexpr);
constexpr double dd_sqr2_sqr_to_double_constexpr = to_double(dd_sqrt2_constexpr*dd_sqrt2_constexpr);

# can't use constinit
double double_stdsqrt2 = sqrt(2);
```

The data section is like

```asm
dd_sqrt2:
        .quad   0x3ff6a09e667f3bcd              # double 1.4142135623730951
        .quad   0xbc9bdd3413b26458              # double -9.6672933134529159E-17

double_sqrt2:
        .quad   0x3ff6a09e667f3bcd              # double 1.4142135623730951

double_sqrt2_sqr:
        .quad   0x4000000000000001              # double 2.0000000000000004

dd_sqrt2_sqr_to_double:
        .quad   0x4000000000000000              # double 2
        
double_stdsqrt2:
        .quad   0x0000000000000000              # double 0
```

Note:

* `fb:function_name` supports compiling time evaluation for `double`
* `qdpp` supports compiling time evaluation for `dd_real`.
* `double_stdsqrt2` will be initialized at the beginning of the execution.
* It's more accuracy  to use `dd_real` as temporary variables.

# Compiler supports

C++20 is need because the use of `std::bit_cast`.

For Clang we use `__buildin_bit_cast` as a workaround.

# Build

* `qdpp` is header only project for c++. 
  * `#include <qd/dd.h>` for all stuff about `dd_real`.
    * For only `dd_real` type definition and basic operations, include `#include <qd/dd_real.h>`.
  * `#include <qd/qd.h>` for all stuff about `qd_real`.
    * For only `qd_real` type definition and basic operations, include `#include <qd/qd_real.h>`.

* If you want to use `qdpp` in C , the simplest way to do that include the `c_*.cpp` in your projects.

  ## Build tests

  `cmake -S ../tests`

# So me notes on changes to `QD`

* `log()` refined for `dd_real`

Refine log implementation for x near 1.

* `read()` refined for `dd_real`

Fix some bugs!
