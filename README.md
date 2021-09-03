

# What is this library?

This library `qdpp` supports double-`double` (104bits) and quad-`double` (209bits) floating points.

This library is a modification of  `QD` library. https://www.davidhbailey.com/dhbsoftware/

# Why is this library?

`qdpp` is based on the `QD` library. 

## So why `QD`?

`QD` supports double-`double` (104bits) and quad-`double` (206bits) floating points. The speed of double-`double` and quad-`double` are much faster than the normal multiple-precision library. Because

* Dynamic memory allocation is not needed. The memory for bits is allocated on stack.
* Easy to implement compiling time evaluation math functions (`constexpr`). 

## `qdpp` compared to `QD`

* Make the library header-only for c++.
* Add `constexpr`  to math functions for `dd_real` and `qd_real`.
* Improvements for some math and IO functions.
* `constexpr` math functions have been implemented for double too.
* `qdpp` is not 100% compatible to the `QD` library.

## Compute constants in compiling time
```c++
constinit dd_real dd_sqrt2 = sqrt(dd_real(2.));

constinit double double_sqrt2 = fb::sqrt(2.);

// 1.4142135623730951
constexpr double double_sqrt2_constexpr = fb::sqrt(2.);
// 2.0000000000000004
constinit double double_sqrt2_sqr = double_sqrt2_constexpr*double_sqrt2_constexpr;

// 1.4142135623730951 -9.6672933134529159E-17
constexpr dd_real dd_sqrt2_constexpr = sqrt(dd_real(2.));;
// 2.
constinit double dd_sqrt2_sqr_to_double = to_double(dd_sqrt2_constexpr*dd_sqrt2_constexpr);

constexpr double dd_sqr2_sqr_to_double_constexpr = to_double(dd_sqrt2_constexpr*dd_sqrt2_constexpr);

// can't use constinit
double double_stdsqrt2 = sqrt(2);
```

You can play with `qdpp` on godbolt  https://godbolt.org/z/xx7eznKTK.

Note:

* ***Functions in namespace `fb` supports compiling time evaluation math functions for `double`***
* ***`qdpp` supports compiling time evaluation for `dd_real/qd_real`.***
* ***It's more accuracy  to use `dd_real` as temporary variables.***

# Limitation

* ***Fast math must be turn off for compiling.***
* ***Floating operation mode must be round-to-nearest.***

# Compiler&Platform supports

GCC & Clang & MSVC are OK to use.

C++20 is needed because of the use of `std::bit_cast`. For Clang we use `__buildin_bit_cast` as a workaround.

Both Windows and Linux (not tested) should be OK.

# Build

* `qdpp` is header only project for c++. 
  
  * All files of c++ parts have been ***clone*** into a single file, `<qd_single.h>`, you can include it for all functions. However, you can also include individual parts only.
    * `#include <qd/double_math.h>` for all stuff about `double`.
    * `#include <qd/dd.h>` for all stuff about `dd_real`.
      * For only `dd_real` type definition and basic operations, include `#include <qd/dd_real.h>`.
    * `#include <qd/dd_real.h>` for all stuff about `qd_real`. 
  
* If you want to use `qdpp` in C , the simplest way to do that include the source files`c_*.cpp` in your projects.

  ## Build tests (Cmake)

  ```bash
  mkdir build
  cd build
  cmake -S ../tests
  ```

# Reference

(Sorry!  Far away to complete.)

#### `double`

Any function in namespace `fb` is implemented like this

```c++
double foo() noexcept {
    if(std::is_constant_evaluated()) {
        my_foo();
    } else {
        return std::foo();
    }
}
```

It doesn't touch `errno` or throw exceptions in constant evaluation context. If any error occur, It just returns `nan`.

### Precision

See the precision.txt for details.

#### `double`

* `round`, `ceil`, `floor`, `abs`, and `fmod` should be exact.
* `sqrt` should exactly round to nearest. (ups <= 1/2 ups)
* `exp`, `log`, `log10`, `pow`, `sin`, `cos` , `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh` has an error up to (0.5+tiny number) ups. Almost all results round correctly to nearest.

#### `qd_real` and `dd_real`

For  basic operations like `+`,`-`,`*`, `/`, the errors are about 0.5-2 ups for qd_real and 1-4 ups for dd_ral.

For most math functions (except `fmod`, `sin`, `cos` and `tan`), rounding errors are < 10 ups for qd_real and <20 ups for dd_real.

However, the accuracy for `fmod`, `sin`, `cos`, and `tan` can be very bad. The rounding error can be thousands  or millions ups.

# Some notes on changes to `QD`

(Not completed)

* Improve the accuracy for various math functions.
* Add `expm1` and `log1p`.

* `ddrand` and `qdrand` now need parameter of generator, `lrand48` is not portable, and `std::rand()` is too bad. Add `real_rand` function template for `ddrand` and `qdrand`
* Remove `*d_real::*_real(char const *)`. Add string literals `_dd` and `_qd`.
* Test the accuracy.
* Test the performance for more operations. `Use `<chrono>` for timing. It's portable, precision, and performant.
