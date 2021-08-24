/*  dd.h  */
#ifndef _QD_DD_H
#define _QD_DD_H

/*  dd_real.h  */
/*
 * include/dd_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Double-double precision (>= 106-bit significand) floating point
 * arithmetic package based on David Bailey's Fortran-90 double-double
 * package, with some changes. See  
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *   
 * for the original Fortran-90 version.
 *
 * Overall structure is similar to that of Keith Brigg's C++ double-double
 * package.  See  
 *
 *   http://www-epidem.plansci.cam.ac.uk/~kbriggs/doubledouble.html
 *
 * for more details.  In particular, the fix for x86 computers is borrowed
 * from his code.
 *
 * Yozo Hida
 */

#ifndef _QD_DD_REAL_H
#define _QD_DD_REAL_H

#define _CRT_SECURE_NO_WARNINGS

#include <cmath>
#include <iostream>
#include <string>
#include <limits>
/*  qd_config.h  */
/* include/qd/qd_config.h.  Generated from qd_config.h.in by configure.  */
#ifndef _QD_QD_CONFIG_H
#define _QD_QD_CONFIG_H  1

/*  config.h  */
#ifndef QD_CONFIGH_H
#define QD_CONFIGH_H
/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to alternate name for `main' routine that is called from a `main' in
   the Fortran libraries. */
#define FC_MAIN main

/* Define to 1 if your system has the clock_gettime function. */
#define HAVE_CLOCK_GETTIME 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if Fortran interface is to be compiled. */
#define HAVE_FORTRAN 1

/* Define to 1 if you have the <fpu_control.h> header file. */
/* #undef HAVE_FPU_CONTROL_H */

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <ieeefp.h> header file. */
#define HAVE_IEEEFP_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if stdbool.h conforms to C99. */
#define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if the system has the type `_Bool'. */
/* #undef HAVE__BOOL */

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#define LT_OBJDIR ".libs/"

/* qd major version number */
#define MAJOR_VERSION 2

/* qd minor version number */
#define MINOR_VERSION 3

/* Name of package */
#define PACKAGE "qd"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "yozo@cs.berkeley.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "qd"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "qd 2.3.12"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "qd"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.3.12"

/* qd patch number (sub minor version) */
#define PATCH_VERSION 12

/* Any special symbols needed for exporting APIs. */
#define QD_API /**/

/* Define this macro to be the copysign(x, y) function. */
#define QD_COPYSIGN(x, y) ::copysign(x, y)

/* Define to 1 to enable debugging code. */
/* #undef QD_DEBUG */

/* If fused multiply-add is available, define correct macro for using it. */
/* #undef QD_FMA */

/* If fused multiply-subtract is available, define correct macro for using it.
   */
/* #undef QD_FMS */

/* Define to 1 if your compiler have the C++ standard include files. */
#define QD_HAVE_STD 1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Define to 1 if your <sys/time.h> declares `struct tm'. */
/* #undef TM_IN_SYS_TIME */

/* Version number of package */
#define VERSION "2.3.12"

/* Whether to use x86 fpu fix. */
#define X86 1

#endif

/* If fused multiply-add is available, define to correct macro for
   using it.  It is invoked as QD_FMA(a, b, c) to compute fl(a * b + c). 
   If correctly rounded multiply-add is not available (or if unsure), 
   keep it undefined.*/
#ifndef QD_FMA
/* #undef QD_FMA */
#endif

/* If fused multiply-subtract is available, define to correct macro for
   using it.  It is invoked as QD_FMS(a, b, c) to compute fl(a * b - c). 
   If correctly rounded multiply-subtract is not available (or if unsure), 
   keep it undefined.*/
#ifndef QD_FMS
/* #undef QD_FMS */
#endif

/* Set the following to 1 to make the addition and subtraction
   to satisfy the IEEE-style error bound 

      fl(a + b) = (1 + d) * (a + b)

   where |d| <= eps.  If set to 0, the addition and subtraction
   will satisfy the weaker Cray-style error bound

      fl(a + b) = (1 + d1) * a + (1 + d2) * b

   where |d1| <= eps and |d2| eps.           */
#ifndef QD_IEEE_ADD
#define QD_IEEE_ADD 1
#endif

/* Set the following to 1 to use slightly inaccurate but faster
   version of multiplication. */
// #define QD_SLOPPY_MUL 1

/* Set the following to 1 to use slightly inaccurate but faster
   version of division. */
// #define QD_SLOPPY_DIV 1

/* Define this macro to be the isfinite(x) function. */
#ifndef QD_ISFINITE
#define QD_ISFINITE(x) fb::isfinite(x)
#endif

/* Define this macro to be the isinf(x) function. */
#ifndef QD_ISINF
#define QD_ISINF(x) fb::isinf(x)
#endif

/* Define this macro to be the isnan(x) function. */
#ifndef QD_ISNAN
#define QD_ISNAN(x) fb::isnan(x)

#endif

/* For C++ only, inline constexpr */
#define QD_CONSTEXPR constexpr
#define QD_HAS_CONSTEXPR 1

#endif /* _QD_QD_CONFIG_H */

struct dd_real {
  double x[2];

  constexpr dd_real(double hi, double lo) : x() { x[0] = hi; x[1] = lo; }
  constexpr dd_real() : x() {x[0] = 0.0; x[1] = 0.0; }
  constexpr dd_real(double h) : x() { x[0] = h; x[1] = 0.0; }
  constexpr dd_real(int h) : x() {
    x[0] = (static_cast<double>(h));
    x[1] = 0.0;
  }

  constexpr explicit dd_real (const double *d) : x() {
    x[0] = d[0]; x[1] = d[1];
  }

  constexpr double _hi() const { return x[0]; }
  constexpr double _lo() const { return x[1]; }

  static const dd_real _2pi;
  static const dd_real _pi;
  static const dd_real _3pi4;
  static const dd_real _pi2;
  static const dd_real _pi4;
  static const dd_real _e;
  static const dd_real _log2;
  static const dd_real _log10;
  static const dd_real _nan;
  static const dd_real _inf;

  static const double _eps;
  static const double _min_normalized;
  static const dd_real _max;
  static const dd_real _safe_max;
  static const int _ndigits;

  QD_CONSTEXPR bool isnan() const;
  QD_CONSTEXPR bool isfinite() const;
  QD_CONSTEXPR bool isinf() const;

  QD_CONSTEXPR static dd_real add(double a, double b);
  QD_CONSTEXPR static dd_real ieee_add(const dd_real &a, const dd_real &b);
  QD_CONSTEXPR static dd_real sloppy_add(const dd_real &a, const dd_real &b);

  QD_CONSTEXPR dd_real &operator+=(double a);
  QD_CONSTEXPR dd_real &operator+=(const dd_real &a);


  QD_CONSTEXPR static dd_real sub(double a, double b);

  QD_CONSTEXPR dd_real &operator-=(double a);
  QD_CONSTEXPR dd_real &operator-=(const dd_real &a);

  QD_CONSTEXPR dd_real operator+() const;
  QD_CONSTEXPR dd_real operator-() const;

  QD_CONSTEXPR static dd_real mul(double a, double b);

  QD_CONSTEXPR dd_real &operator*=(double a);
  QD_CONSTEXPR dd_real &operator*=(const dd_real &a);

  QD_CONSTEXPR static dd_real div(double a, double b);
  QD_CONSTEXPR static dd_real sloppy_div(const dd_real &a, const dd_real &b);
  QD_CONSTEXPR static dd_real accurate_div(const dd_real &a, const dd_real &b);
  
  QD_CONSTEXPR dd_real &operator/=(double a);
  QD_CONSTEXPR dd_real &operator/=(const dd_real &a);

  QD_CONSTEXPR dd_real &operator=(double a);

  QD_CONSTEXPR dd_real operator^(int n);
  QD_CONSTEXPR static dd_real sqr(double d);

  QD_CONSTEXPR static dd_real sqrt(double a);
  
  QD_CONSTEXPR bool is_zero() const;
  QD_CONSTEXPR bool is_one() const;
  QD_CONSTEXPR bool is_positive() const;
  QD_CONSTEXPR bool is_negative() const;

  static dd_real rand(void);

  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int len, int precision = _ndigits, 
      bool showpos = false, bool uppercase = false) const;
  std::string to_string(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;

  static constexpr int read(const char *s, dd_real &a);
  static constexpr dd_real read(char const* s);

  /* Debugging Methods */
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "", 
                 std::ostream &os = std::cerr) const;

  static dd_real debug_rand();
};
/*  dd_const.inl.h  */
/*
 * src/dd_const.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 */

/*  two_basics.h  */
/*
 * include/inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This file contains the basic functions used both by double-double
 * and quad-double package.  These are declared as inline functions as
 * they are the smallest building blocks of the double-double and 
 * quad-double arithmetic.
 */
#ifndef _QD_INLINE_H
#define _QD_INLINE_H

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

#include <cmath>
#include <limits>
#include <type_traits>

namespace qd {

    inline QD_CONSTEXPR const double _d_nan = std::numeric_limits<double>::quiet_NaN();
    inline QD_CONSTEXPR const double _d_inf = std::numeric_limits<double>::infinity();


#ifdef __GNUC__
#if defined(__FAST_MATH__) && __FAST_MATH__

#include <stdio.h>
#include <type_traits>
#include <exception>
#define TWO_CHECK() do {\
    if (!std::is_constant_evaluated()) {\
        fprintf(stderr, "Don't use fast-math when you use dd_real/qd_real\n");\
        std::terminate();\
    } else {\
        int dont_use_fast_math_when_you_use_dd_real = 1;\
        dont_use_fast_math_when_you_use_dd_real -= dont_use_fast_math_when_you_use_dd_real;\
        int dumy = 1 / dont_use_fast_math_when_you_use_dd_real;\
    }\
} while(false)
#else
#define TWO_CHECK() do {} while(false)
#endif
#else
#define TWO_CHECK() do {} while(false)
#endif


#ifdef _MSC_VER

#include <immintrin.h>

    inline double QD_FMA_NOCHECK(double a, double b, double c)
    {
        __m128d aw;
        aw.m128d_f64[0] = a;
        aw.m128d_f64[1] = 0;
        __m128d bw;
        bw.m128d_f64[0] = b;
        bw.m128d_f64[1] = 0;
        __m128d cw;
        cw.m128d_f64[0] = c;
        cw.m128d_f64[1] = 0;
        __m128d answ = _mm_fmadd_sd(aw, bw, cw);
        return answ.m128d_f64[0];
    }

    inline double QD_FMS_NOCHECK(double a, double b, double c)
    {
        __m128d aw;
        aw.m128d_f64[0] = a;
        aw.m128d_f64[1] = 0;
        __m128d bw;
        bw.m128d_f64[0] = b;
        bw.m128d_f64[1] = 0;
        __m128d cw;
        cw.m128d_f64[0] = c;
        cw.m128d_f64[1] = 0;
        __m128d answ = _mm_fmsub_sd(aw, bw, cw);
        return answ.m128d_f64[0];
    }

#ifndef QD_FMA
#define QD_FMA QD_FMA_NOCHECK
#endif

#ifndef QD_FMS
#define QD_FMS QD_FMS_NOCHECK
#endif

#endif


    /*********** Basic Functions ************/
    /* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. or |a||b| = 0. */
    inline QD_CONSTEXPR double quick_two_sum(double a, double b, double& err)
    {
        TWO_CHECK();
        double s = a + b;
        err = b - (s - a);
        return s;
    }

    /* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| or |a||b| = 0. */
    inline QD_CONSTEXPR double quick_two_diff(double a, double b, double& err)
    {
        TWO_CHECK();
        double s = a - b;
        err = (a - s) - b;
        return s;
    }

    /* Computes fl(a+b) and err(a+b).  */
    inline QD_CONSTEXPR double two_sum(double a, double b, double& err)
    {
        TWO_CHECK();
        double s = a + b;
        double bb = s - a;
        err = (a - (s - bb)) + (b - bb);
        return s;
    }

    /* Computes fl(a-b) and err(a-b).  */
    inline QD_CONSTEXPR double two_diff(double a, double b, double& err)
    {
        TWO_CHECK();
        double s = a - b;
        double bb = s - a;
        err = (a - (s - bb)) - (b + bb);
        return s;
    }

    /* Computes high word and lo word of a */
    inline QD_CONSTEXPR  void split(double a, double& hi, double& lo)
    {
        double temp;
        if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
            a *= 3.7252902984619140625e-09;  // 2^-28
            temp = _QD_SPLITTER * a;
            hi = temp - (temp - a);
            lo = a - hi;
            hi *= 268435456.0;          // 2^28
            lo *= 268435456.0;          // 2^28
        } else {
            temp = _QD_SPLITTER * a;
            hi = temp - (temp - a);
            lo = a - hi;
        }
    }

    /* Computes fl(a*b) and err(a*b). */
    inline QD_CONSTEXPR double two_prod(double a, double b, double& err)
    {
        TWO_CHECK();

#ifdef QD_FMS
        bool const has_fms = true;
#else
        bool const has_fms = false;
#endif

        if (!std::is_constant_evaluated() && has_fms) {
            double p = a * b;
            err = QD_FMS(a, b, p);
            return p;
        } else {
            double a_hi, a_lo, b_hi, b_lo;
            double p = a * b;
            split(a, a_hi, a_lo);
            split(b, b_hi, b_lo);
            err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
            return p;
        }
    }

    /* Computes fl(a*a) and err(a*a).  Faster than the above method. */
    inline QD_CONSTEXPR double two_sqr(double a, double& err)
    {
        TWO_CHECK();

#ifdef QD_FMS
        bool const has_fms = true;
#else
        bool const has_fms = false;
#endif

        if (!std::is_constant_evaluated() && has_fms) {
            double p = a * a;
            err = QD_FMS(a, a, p);
            return p;
        } else {
            double hi, lo;
            double q = a * a;
            split(a, hi, lo);
            err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
            return q;
        }

    }
}

#endif /* _QD_INLINE_H */
inline constexpr dd_real dd_real::_2pi = dd_real(6.283185307179586232e+00,
                                      2.449293598294706414e-16);
inline constexpr  dd_real dd_real::_pi = dd_real(3.141592653589793116e+00,
                                     1.224646799147353207e-16);
inline constexpr dd_real dd_real::_pi2 = dd_real(1.570796326794896558e+00,
                                      6.123233995736766036e-17);
inline constexpr dd_real dd_real::_pi4 = dd_real(7.853981633974482790e-01,
                                      3.061616997868383018e-17);
inline constexpr dd_real dd_real::_3pi4 = dd_real(2.356194490192344837e+00,
                                       9.1848509936051484375e-17);
inline constexpr dd_real dd_real::_e = dd_real(2.718281828459045091e+00,
                                    1.445646891729250158e-16);
inline constexpr dd_real dd_real::_log2 = dd_real(6.931471805599452862e-01,
                                       2.319046813846299558e-17);
inline constexpr dd_real dd_real::_log10 = dd_real(2.302585092994045901e+00,
                                        -2.170756223382249351e-16);
inline constexpr dd_real dd_real::_nan = dd_real(qd::_d_nan, qd::_d_nan);
inline constexpr dd_real dd_real::_inf = dd_real(qd::_d_inf, qd::_d_inf);

inline constexpr double dd_real::_eps = 4.93038065763132e-32;  // 2^-104
inline constexpr double dd_real::_min_normalized = 2.0041683600089728e-292;  // = 2^(-1022 + 53)
inline constexpr dd_real dd_real::_max =
    dd_real(1.79769313486231570815e+308, 9.97920154767359795037e+291);
inline constexpr dd_real dd_real::_safe_max =
    dd_real(1.7976931080746007281e+308, 9.97920154767359795037e+291);
inline constexpr int dd_real::_ndigits = 31;



namespace std {
    template <>
    class numeric_limits<dd_real> : public numeric_limits<double> {
    public:
        static constexpr double epsilon() { return dd_real::_eps; }
        static constexpr dd_real(max)() { return dd_real::_max; }
        static constexpr dd_real safe_max() { return dd_real::_safe_max; }
        static constexpr double (min)() { return dd_real::_min_normalized; }
        static constexpr const int digits = 104;
        static constexpr const int digits10 = 31;
    };
}

QD_API constexpr dd_real operator+(const dd_real &a, double b);
QD_API constexpr dd_real operator+(double a, const dd_real &b);
QD_API constexpr dd_real operator+(const dd_real &a, const dd_real &b);

QD_API constexpr dd_real operator-(const dd_real &a, double b);
QD_API constexpr dd_real operator-(double a, const dd_real &b);
QD_API constexpr dd_real operator-(const dd_real &a, const dd_real &b);

QD_API constexpr dd_real operator*(const dd_real &a, double b);
QD_API constexpr dd_real operator*(double a, const dd_real &b);
QD_API constexpr dd_real operator*(const dd_real &a, const dd_real &b);

QD_API constexpr dd_real operator/(const dd_real &a, double b);
QD_API constexpr dd_real operator/(double a, const dd_real &b);
QD_API constexpr dd_real operator/(const dd_real &a, const dd_real &b);

QD_API constexpr dd_real inv(const dd_real &a);

constexpr bool operator==(const dd_real &a, double b);
constexpr bool operator==(double a, const dd_real &b);
constexpr bool operator==(const dd_real &a, const dd_real &b);

QD_API constexpr bool operator<=(const dd_real &a, double b);
QD_API constexpr bool operator<=(double a, const dd_real &b);
QD_API constexpr bool operator<=(const dd_real &a, const dd_real &b);

QD_API constexpr bool operator>=(const dd_real &a, double b);
QD_API constexpr bool operator>=(double a, const dd_real &b);
QD_API constexpr bool operator>=(const dd_real &a, const dd_real &b);

QD_API constexpr bool operator<(const dd_real &a, double b);
QD_API constexpr bool operator<(double a, const dd_real &b);
QD_API constexpr bool operator<(const dd_real &a, const dd_real &b);

QD_API constexpr bool operator>(const dd_real &a, double b);
QD_API constexpr bool operator>(double a, const dd_real &b);
QD_API constexpr bool operator>(const dd_real &a, const dd_real &b);

QD_API constexpr bool operator!=(const dd_real &a, double b);
QD_API constexpr bool operator!=(double a, const dd_real &b);
QD_API constexpr bool operator!=(const dd_real &a, const dd_real &b);


/* Round to double */
constexpr double to_double(const dd_real& a);
constexpr int    to_int(const dd_real& a);

constexpr dd_real sqr(const dd_real& a);
constexpr dd_real npwr(const dd_real& a, int n);
/* Computes  dd * d  where d is known to be a power of 2. */
constexpr dd_real mul_pwr2(const dd_real& dd, double d);



/*  dd_inline.h  */
/*
 * include/dd_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains small functions (suitable for inlining) in the double-double
 * arithmetic package.
 */
#ifndef _QD_DD_INLINE_H
#define _QD_DD_INLINE_H

#include <cmath>
/*  double_basics.h  */
#ifndef QD_FLOATING_BASICS_H
#define QD_FLOATING_BASICS_H

#include <stdint.h>
#include <utility>
#include <bit>
#include <cmath>
#include <limits>
#include <assert.h>

namespace fb {

    //         double(_pi)      = 3.14159265358979311599796346854418516159057617187500
    inline constexpr double _pi = 3.1415926535897932384626433832795028841971693993751058209749445923;
    inline constexpr double _d_pi_t1 = 1.2246467991473531772260659322750011E-16;
    // pi ~ _pi + _pi_t1
    //                                   1.57079632679489655799898173427209258079528808593750
    inline constexpr double _d_half_pi = 1.5707963267948966192313216916397514420985846996875529104874722961;
    inline constexpr double _d_pi = 3.1415926535897932384626433832795028841971693993751058209749445923;
    inline constexpr double _d_2pi = 6.2831853071795864769252867665590057683943387987502116419498891846;

    inline constexpr double _eps = std::numeric_limits<double>::epsilon();
    inline constexpr double _d_nan = std::numeric_limits<double>::quiet_NaN();
    inline constexpr double _subnorm_min = 4.9406564584124654E-324;
    inline constexpr double _subnorm_max = 2.2250738585072009E-308;
    inline constexpr double _norm_min = 2.2250738585072014E-308;
    inline constexpr double _max = 1.7976931348623157E+308;
    inline constexpr double _infinity = INFINITY;
    inline constexpr double _NaN = NAN;
    inline constexpr double _int_min = -9007199254740992;
    inline constexpr double _int_max = +9007199254740992;
    inline constexpr double _d_zero = 0.0;
    inline constexpr double _d_one = 1.0;
    inline constexpr double _d_one_next = 1.0000000000000002;
    inline constexpr double _d_one_prev = 0.9999999999999999;
    inline constexpr double _d_two = 2.0;
    //                               0.69314718055994528622676398299518041312694549560546875
    inline constexpr double _d_ln2 = 0.6931471805599453094172321214581765680755001343602552541206800094;
    inline constexpr double _d_ln2_t1 = 2.319046813846299615494855463875478650E-17;
    inline constexpr double _d_ln3 = 1.0986122886681096913952452369225257046474905578227494517346943336;
    inline constexpr double _d_ln4 = 1.3862943611198906188344642429163531361510002687205105082413600189;
    inline constexpr double _d_ln5 = 1.6094379124341003746007593332261876395256013542685177219126478914;
    inline constexpr double _d_ln6 = 1.7917594692280550008124773583807022727229906921830047058553743431;
    inline constexpr double _d_ln7 = 1.9459101490553133051053527434431797296370847295818611884593901499;
    inline constexpr double _d_ln10 = 2.3025850929940456840179914546843642076011014886287729760333279009;

    inline constexpr double _d_ln_pi =   1.1447298858494001741434273513530587116472948129153115715136230714;
    inline constexpr double _d_ln_8pi = 3.2241714275292361023951237157275884158737952159960773338756630999;

    inline constexpr double _d_sin_1 = 0.8414709848078965066525023216302989996225630607983710656727517099;
    inline constexpr double _d_sin_1d5 = 0.9974949866040544309417233711414873227066514259221158219499748240;
    inline constexpr double _d_e = 2.7182818284590452353602874713526624977572470936999595749669676277;
    inline constexpr uint64_t to_uint64(double x)
    {
#ifdef _MSC_VER
        return std::bit_cast<uint64_t>(x);
#else
        return __builtin_bit_cast(uint64_t, x);
#endif
    }

    inline constexpr double to_double(uint64_t x)
    {
#ifdef _MSC_VER
        return std::bit_cast<double>(x);
#else
        return __builtin_bit_cast(double, x);
#endif
    }

    constexpr int EXP_RAW_INFINITY = 0x7FF;
    constexpr int EXP_RAW_ERROR = 0x7FF;

    constexpr int EXP_MAX = 1023;
    constexpr int EXP_NORMAL_MIN = -1022;
    constexpr int EXP_SUBNORMAL = -1023;

    constexpr int EXP_MIN = -1074;
    constexpr int EXP_OFFSET = 52;

    constexpr uint64_t FRAC_MASK = (uint64_t(1) << EXP_OFFSET) - 1;
    constexpr int EXP_BIAS = 1023;
    constexpr int SIGN_OFFSET = 63;


    struct Float128 {
    };

    constexpr uint64_t full_number(uint64_t x)
    {
        return x | (uint64_t(1) << EXP_OFFSET);
    }

    inline constexpr int exp_part(double u)
    {
        uint64_t u64 = to_uint64(u);
        return int((u64 >> EXP_OFFSET) & 0x7FF) - EXP_BIAS;
    }

    inline constexpr int raw_exp(double u)
    {
        return (to_uint64(u) >> EXP_OFFSET) & 0x7FF;
    }

    inline constexpr uint64_t frac_part(double u)
    {
        return to_uint64(u) & FRAC_MASK;
    }

    inline constexpr uint64_t frac_part(uint64_t u)
    {
        return u & FRAC_MASK;
    }

    inline constexpr double make_double(uint64_t frac, int exp, bool sign)
    {
        uint64_t v = frac | (uint64_t(exp + EXP_BIAS) << EXP_OFFSET)
            | (uint64_t(sign) << SIGN_OFFSET);
        return to_double(v);
    }


    // zero means the least significant bit
    inline constexpr int leading_bit(uint64_t u)
    {
        for (int i = 63; i >= 0; --i) {
            if (u & (uint64_t(1) << i)) {
                return i;
            }
        }
        return -1;

    }
    static_assert(leading_bit(0) == -1);
    static_assert(leading_bit(1) == 0);
    static_assert(leading_bit(-1) == 63);

    inline constexpr bool isinf_(double x)
    {
        uint64_t u = to_uint64(x);
        return raw_exp(x) == EXP_RAW_INFINITY && frac_part(x) == 0;
    }

    inline constexpr bool isinf(double x)
    {
        if (std::is_constant_evaluated()) {
            return isinf_(x);
        } else {
            return std::isinf(x);
        }
    }

    inline constexpr bool isnan_(double x)
    {
        uint64_t u = to_uint64(x);
        return raw_exp(x) == EXP_RAW_INFINITY && frac_part(x) != 0;
    }

    inline constexpr bool isnan(double x)
    {
        if (std::is_constant_evaluated()) {
            return isnan_(x);
        } else {
            return std::isnan(x);
        }
    }

    inline constexpr bool isfinite_(double x)
    {
        return raw_exp(x) != EXP_RAW_INFINITY;
    }

    inline constexpr bool isfinite(double x)
    {
        if (std::is_constant_evaluated()) {
            return isfinite_(x);
        } else {
            return std::isfinite(x);
        }
    }

    inline constexpr bool signbit_(double x)
    {
        return to_uint64(x) & (uint64_t(1) << 63);
    }

    inline constexpr bool signbit(double x)
    {
        if (std::is_constant_evaluated()) {
            return signbit_(x);
        } else {
            return std::signbit(x);
        }
    }

    inline constexpr double copysign_(double x, double sign)
    {
        uint64_t bit = (uint64_t)signbit(sign) << 63;
        uint64_t const mask = (uint64_t(1) << 63) - uint64_t(1);
        return to_double((to_uint64(x) & mask) | bit);
    }

    inline constexpr double copysign(double x, double sign)
    {
        if (std::is_constant_evaluated()) {
            return copysign_(x, sign);
        } else {
            return std::copysign(x, sign);
        }
    }

    inline constexpr double abs_(double x)
    {
        if (std::is_constant_evaluated()) {
            return copysign(x, 0);
        } else {
            return std::abs(x);
        }
    }

    inline constexpr double abs(double x)
    {
        if (std::is_constant_evaluated()) {
            return abs_(x);
        } else {
            return std::abs(x);
        }
    }
    
    // if x == INT_MAX return -x
    inline constexpr int abs_(int x)
    {
        if (x < 0) return -x;
        else return x;
    }

    inline constexpr int abs(int x)
    {
        if (std::is_constant_evaluated()) {
            return abs_(x);
        } else {
            return std::abs(x);
        }
    }

    inline constexpr double sqr_(double d)
    {
        return d * d;
    }

    inline constexpr double floor__(double d)
    {
        if (isnan_(d) || isinf_(d) || d == 0.) {
            // if d == +=0. || d is inf || d is nan
            //  return d
            return d;
        }
        bool sign = d < 0;
        d = abs_(d);
        if (sign && d < 1.) {
            return -1;
        } else if (!sign && d < 1.) {
            return 0.;
        } else if (d >= 4503599627370496) { // 2^52
            return sign ? -d : d;
        } else {
            uint64_t di = uint64_t(d);
            if (di == d) {
                return sign ? -d : d;
            } else {
                if (!sign) return (double)di;
                else return -(double)(di + 1);
            }
        }
    }

    inline constexpr double floor_(double d)
    {
        if (isnan_(d) || isinf_(d) || d == 0.) {
            // if d == +=0. || d is inf || d is nan
            //  return d
            return d;
        }
        bool sign = signbit(d) ? true : false;
        d = abs(d);

        int exp_part_ = exp_part(d);
        if (sign && d < 1.)
            return -1.;
        else if (!sign && d < 1.) {
            return 0.;
        } else if (exp_part_ >= 52) { // already integer
            return sign ? -d : d;
        } else { // normal integer
            uint64_t double_frac_part = frac_part(d);
            uint64_t const double_frac_mask = (uint64_t(1) << 52) - 1;
            uint64_t int_bit = uint64_t(1) << (52 - exp_part_);
            uint64_t int_mask = ~(int_bit - 1);
            uint64_t frac_mask = int_bit - 1;
            uint64_t frac_part_2 = double_frac_part & int_mask;
            if (sign) {
                // add leading one
                frac_part_2 = frac_part_2 | (uint64_t(1) << 52);

                // ceil
                if (double_frac_part & frac_mask) {
                    frac_part_2 += int_bit;
                }
                if (frac_part_2 >= (uint64_t(1) << 53)) {
                    exp_part_ += 1;
                    frac_part_2 >>= 1;
                }

                // erase leading one
                frac_part_2 &= double_frac_mask;
            }
            double v = make_double(frac_part_2, exp_part_, sign);
            return v;
        }
    }

    inline constexpr double floor(double d)
    {
        if (std::is_constant_evaluated()) {
            return floor_(d);
        } else {
            return std::floor(d);
        }
    }

    inline constexpr double ceil_(double d)
    {
        if (isnan(d) || isinf(d) || d == 0.) {
            // if d == +=0. || d is inf || d is nan
            //  return d
            return d;
        }
        // finite
        return -floor_(-d);
    }

    inline constexpr double ceil(double d)
    {
        if (std::is_constant_evaluated()) {
            return ceil_(d);
        } else {
            return std::ceil(d);
        }
    }

    inline constexpr double round_(double d)
    {
        if (isnan_(d) || isinf_(d) || d == 0.) {
            return d;
        }
        bool sign = signbit_(d);
        d = abs_(d);
        if (d < 0.5)
            return sign ? -0. : 0.;
        else {
            return sign ? -floor_(d + 0.5) : floor_(d + 0.5);
        }
    }

    inline constexpr double round(double d)
    {
        if (std::is_constant_evaluated()) {
            return round_(d);
        } else {
            return std::round(d);
        }
    }

    // trunc the result if submnorm
    inline constexpr double ldexp_(double x, int exp)
    {
        if (x == 0. || isinf(x) || isnan(x) || exp == 0) {
            return x;
        } else {
            uint64_t frac = frac_part(x);
            int exp_part_ = exp_part(x);

            if (exp_part_ >= EXP_NORMAL_MIN) {
                // normal               
                exp_part_ += exp;
            } else {
                // subnormal
                int leading_bit_ = leading_bit(frac);
                frac <<= (52 - leading_bit_);
                frac &= FRAC_MASK;
                exp_part_ = EXP_SUBNORMAL + 1 + exp - (52 - leading_bit_);
            }


            // we first check if there is an overflow for adding in (exp + exp_part_)
            if (exp > 2048 || exp_part_ >= 1024) { // overflow
                return copysign(INFINITY, x);
            } else if (exp < -2048 || exp_part_ < -1074) {
                return copysign(0., x);
            } else if (exp_part_ >= EXP_NORMAL_MIN) { // normal
                double v = make_double(frac, exp_part_, signbit(x));
                return v;
            } else { // subnormal
                int shift = (1 + EXP_SUBNORMAL - exp_part_);
                frac |= uint64_t(1) << EXP_OFFSET;
                frac >>= shift;
                double v = make_double(frac, EXP_SUBNORMAL, signbit(x));
                return v;
            }
        }
    }

    inline constexpr double ldexp(double x, int exp)
    {
        if (std::is_constant_evaluated()) {
            return ldexp_(x, exp);
        } else {
            return std::ldexp(x, exp);
        }
    }

    inline constexpr double frexp_(double arg, int* exp)
    {
        if (arg == 0.) {
            *exp = 0;
            return arg;
        } else if (isinf(arg)) {
            return arg;
        } else if (isnan(arg)) {
            return arg;
        } else {
            int exp_ = exp_part(arg);
            if (exp_ != EXP_SUBNORMAL) {
                exp_ += 1;
            } else { // subnormal
                uint64_t u = to_uint64(arg);
                int leading_bit_ = leading_bit(u);
                exp_ = -1022 - (51 - leading_bit_);
            }
            *exp = exp_;
            return ldexp_(arg, -exp_);
        }
    }

    inline constexpr double frexp(double arg, int* exp)
    {
        if (std::is_constant_evaluated()) {
            return frexp_(arg, exp);
        } else {
            return std::frexp(arg, exp);
        }
    }
    // round the result if submnorm
    inline constexpr double mul2pwr_(double x, int exp)
    {
        if (x == 0. || isinf(x) || isnan(x) || exp == 0) {
            return x;
        } else {
            int a;
            double x_ = frexp_(x, &a);
            x_ *= 2;
            a -= 1;

            a += exp;
            if (a >= EXP_MIN) {
                double v = ldexp(1., a);
                x_ *= v;
                return x_;
            } else if (a == EXP_MIN - 1) {
                return _subnorm_min;
            } else {
                return copysign(0, x);
            }
        }
    }

    struct Frexp2 { double x; int exp; };
    inline constexpr Frexp2 frexp_2(double arg)
    {
        Frexp2 ans = { };
        if (arg == 0.) {
            ans.x = arg;
        } else if (isinf_(arg)) {
            ans.x = arg;
        } else if (isnan_(arg)) {
            ans.x = arg;
        } else {
            int exp_ = exp_part(arg);
            if (exp_ != EXP_SUBNORMAL) {
                exp_ += 1;
            } else { // subnormal
                uint64_t u = to_uint64(arg);
                int leading_bit_ = leading_bit(u);
                exp_ = -1022 - (51 - leading_bit_);
            }
            ans.exp = exp_;
            ans.x = ldexp_(arg, -exp_);
        }
        return ans;
    }

    struct Uint128 {
        uint64_t m_hig;
        uint64_t m_low;

        constexpr Uint128(uint64_t low)
            : m_low(low), m_hig(0)
        {
        }

        constexpr Uint128() : Uint128(0, 0)
        {
        }

        constexpr Uint128(uint64_t hig, uint64_t low)
            : m_hig(hig), m_low(low)
        {
        }

        static constexpr Uint128 mul(uint64_t u, uint64_t v)
        {
            constexpr const uint32_t uint32_mask = -1;
            Uint128 ans(0);
            ans.m_low += (u & uint32_mask) * (v & uint32_mask);
            ans.m_hig += ((u>>32) & uint32_mask) * ((v>>32) & uint32_mask);

            uint64_t m1 = ((u >> 32) & uint32_mask) * (v & uint32_mask);
            uint64_t m2 = (u & uint32_mask) * ((v >> 32) & uint32_mask);

            ans = Uint128::add(ans, (m1 & uint32_mask) << 32);
            ans = Uint128::add(ans, (m2 & uint32_mask) << 32);

            ans.m_hig += (m1 >> 32);
            ans.m_hig += (m2 >> 32);
            return ans;
        }

        static constexpr Uint128 mul(Uint128 u, uint64_t v)
        {
            Uint128 p1 = Uint128::mul(u.m_low, v);
            Uint128 p2 = Uint128::mul(u.m_hig, v);
            p1.m_hig += p2.m_low;
            return p1;
        }

        static constexpr Uint128 mul(Uint128 u, Uint128 v)
        {
            Uint128 p1 = Uint128::mul(u, v.m_low);
            Uint128 p2 = Uint128::mul(u, v.m_hig);
            p1.m_hig += p2.m_low;
            return p1;
        }

        static constexpr Uint128 add(uint64_t u, uint64_t v)
        {
            return Uint128((u + v) < u, u + v);
        }
        static constexpr Uint128 add(Uint128 u, uint64_t v)
        {
            u.m_low += v;
            if (u.m_low < v) {
                u.m_hig += 1;
            }
            return u;
        }
        static constexpr Uint128 add(Uint128 u, Uint128 v)
        {
            Uint128 p1 = Uint128::add(u, v.m_low);
            p1.m_hig += v.m_hig;
            return p1;
        }
    };

    inline constexpr bool operator<(Uint128 lhs, Uint128 rhs)
    {
        return lhs.m_hig < rhs.m_hig ||
            (lhs.m_hig == rhs.m_hig && lhs.m_low < rhs.m_low);
    }

    inline constexpr bool operator<=(Uint128 lhs, Uint128 rhs)
    {
        return lhs.m_hig < rhs.m_hig ||
            (lhs.m_hig == rhs.m_hig && lhs.m_low <= rhs.m_low);
    }

    inline constexpr bool operator==(Uint128 lhs, Uint128 rhs)
    {
        return lhs.m_hig == rhs.m_hig && lhs.m_low == rhs.m_low;
    }

    static_assert(Uint128::mul(1, 1) == Uint128(0, 1));
    static_assert(Uint128::mul(uint64_t(1) << 32, 1) == Uint128(0, uint64_t(1) << 32));
    static_assert(Uint128::mul(uint64_t(1) << 32, 2) == Uint128(0, uint64_t(1) << 33));
    static_assert(Uint128::mul(1, uint64_t(1) << 32) == Uint128(0, uint64_t(1) << 32));
    static_assert(Uint128::mul(2, uint64_t(1) << 32) == Uint128(0, uint64_t(1) << 33));
    static_assert(Uint128::mul(uint64_t(1) << 32, uint64_t(1) << 32) == Uint128(1, 0));
    static_assert(Uint128::mul(uint64_t(1) << 63, uint64_t(1) << 63) == Uint128(uint64_t(1) << (63 - 1), 0));
    // (2^64-1)^2 = 2^128 - 2 * 2^64 + 1
    static_assert(Uint128::mul(uint64_t(-1), uint64_t(-1)) == Uint128(-2, 1));

    constexpr Uint128 v = Uint128::mul(uint64_t(-1), uint64_t(-1));

    // sqrt((1<<64) * d)
    // fraction part truced
    inline constexpr uint64_t sqrt_shift(uint64_t d)
    {
        Uint128 u128(d, 0);
        uint64_t x = 0;
        for (int i = 63; i >= 0; --i) {
            uint64_t test = (uint64_t(1) << i) | x;
            if (Uint128::mul(test, test) <= u128) {
                x |= (uint64_t(1) << i);
            }
        }
        return x;
    }

    static_assert(sqrt_shift(0) == 0);
    static_assert(sqrt_shift(1) == uint64_t(1) << 32);
    // sqrt(A^2 - 1) ~ A (1-1/A^2) = A - 1/A
    static_assert(sqrt_shift(18014398509481982) == uint64_t(576460752303423455));

    // sqrt((1<<64) * d)
    // fraction part rounded
    inline constexpr uint64_t sqrt_shift_round(uint64_t d)
    {
        uint64_t x = sqrt_shift(d);
        // x*x <= v
        // (x+1)*(x+1) > v
        // (x + 0.5)^2 = x^x + x + 0.25
        Uint128 int_part = Uint128::add(Uint128::mul(x, x), x);
        Uint128 v = Uint128(d, 0);
        if (v <= int_part) {
            return x;
        } else {
            // may overflow
            return x + 1;
        }
    }

    static_assert(sqrt_shift_round(0) == 0);
    static_assert(sqrt_shift_round(1) == uint64_t(1) << 32);
    // (2^64 - 0.5)^2 =  2^128 - 2^64 + 0.25
    // sqrt(2^64 (2^64-1)) = 2^64 - 1
    static_assert(sqrt_shift_round(-1) == uint64_t(-1));



    inline constexpr double sqrt_(double d)
    {
        if (isnan_(d) || d == 0.)
            return d;
        else if (d < 0.) {
            return std::numeric_limits<double>::quiet_NaN();
        } else if (isinf_(d)) {
            return d;
        } else {
            // positive and finite
            int exp_part_ = exp_part(d);
            if (exp_part_ == -1023) { // subnormal
                // 0.frac 2^-1022
                double constexpr pre_mul = double(uint64_t(1) << 55);
                return sqrt(pre_mul * pre_mul * d) / pre_mul;
            } else { // normal
                // 1.frac 2^e
#if 1
                uint64_t frac = frac_part(d);
                bool mul2 = exp_part_ & 1;
                int exp_part_half = exp_part_ >> 1;

                uint64_t number = full_number(frac);  // [1 - 2)
                if (mul2) number *= 2;                // [1 - 4)
                // round then round is not correct
                // e.g. 4.49 - round -> 4.5 round 5.
                uint64_t value = sqrt_shift(number);  // [1 - 2)


                // number = real value << (52 + 64) (ignore the exp)
                // we need shift right (52 + 64) / 2 then shift left 52
                // we do it simutaneously to avoid losing accuracy                
                const int shift_base = ((52 + 64) / 2 - 52);

                if ((value >> shift_base) & (uint64_t(1) << 53)) {
                    assert(false);
                }

                frac = value >> shift_base;

                if ((value >> (shift_base-1)) & unsigned(1)) {
                    frac += 1;
                    if ((uint64_t(1) << 53) & frac) {
                        exp_part_half += 1;
                        frac >>= 1;
                    }
                } // else fraction truc-ed
                frac &= FRAC_MASK;
                double x = make_double(frac, exp_part_half, signbit(d));
                return x;

#else
                    // (2)*1.frac 2^{2*m}
                    // sqrt() = sqrt(2*1.frac) 2^m
                    // sqrt() ~ 1.414 * (1 + 0.frac/2) 2^m
                bool mul2 = exp_part_ & 1;
                int exp_part_half = exp_part_ >> 1;
                uint64_t frac = frac_part(d);
                frac /= 2;

                if (!frac && !mul2) { // pow of four
                    return make_double(frac, exp_part_half, signbit(d));
                }

                double x = make_double(frac, exp_part_half, signbit(d));
                if (mul2) x *= 1.41421;

                // solve x for f(x) = y
                // x = x - (f(x) - y)/f'(x)
                // f'(x) = x^2' = 2x
                // x = x - (x^2 - y)/(2x)
                // x = 0.5*(x + y/x)

                double old_x = -0.;
                for (int i = 0; i < 100 && (old_x != x); ++i) {
                    old_x = x;
                    // round error!!!!!!!!!!!!!!!!!!!!
                    x = 0.5 * (x + d / x);
                }
                return x;
#endif

            }

        }

    }
    inline constexpr double sqrt(double d) noexcept
    {
        if (std::is_constant_evaluated()) {
            return sqrt_(d);
        } else {
            return std::sqrt(d);
        }
    }


    inline constexpr double fmod_(double x, double y) noexcept
    {

        if (isnan_(x) || isnan(y)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        //  x\y         0    finite   inf
        //  +-0         nan  +-0      +-0
        //  finite      nan           x
        //  inf         nan  nan      nan

        if (x == 0.) {
            if(!(y == 0.))
                return x;
            else
                return std::numeric_limits<double>::quiet_NaN();
        } else if (isinf_(x)) {
            return std::numeric_limits<double>::quiet_NaN();
        } else {

            if (y == 0.) {
                return std::numeric_limits<double>::quiet_NaN();
            } else if (isinf_(y)) {
                return x;
            } else {
                // both are non-zero finite
                // https://stackoverflow.com/questions/62785780/how-do-i-implement-a-modulus-operator-for-double-variables-using-frexp/62865640#62865640
                double const absx = abs_(x);
                double const absy = abs_(y); // sign of y has no effect

                if (absx >= absy) {
                    double dividend = absx;
                    int expx;
                    int expy;
                    (void)frexp_(absx, &expx);
                    (void)frexp_(absy, &expy);
                    // expx >= expy
                    double divisor = ldexp(absy, expx - expy);
                    // divisor is still finite
                    // 
                    // absx = 0.1xxxx 2^expx
                    // absy = 0.1yyyy 2^expx
                    // 2*absy > absx

                    for (;divisor >= absy;) {
                        // dividend < 2*divisor
                        if (dividend >= divisor) {
                            dividend -= divisor;
                        }
                        divisor *= 0.5;
                    }
                    return copysign_(dividend, x);
                } else {
                    return x;
                }
            }
        }
    }

    inline constexpr double fmod(double x, double y) noexcept
    {
        if (std::is_constant_evaluated()) {
            return fmod_(x, y);
        } else {
            return std::fmod(x, y);
        }
    }

    inline constexpr double remainder_(double x, double y) noexcept
    {

        if (isnan_(x) || isnan(y)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        //  x\y         0    finite   inf
        //  +-0         nan  +-0      +-0
        //  finite      nan           x
        //  inf         nan  nan      nan

        if (x == 0.) {
            if (!(y == 0.))
                return x;
            else
                return std::numeric_limits<double>::quiet_NaN();
        } else if (isinf_(x)) {
            return std::numeric_limits<double>::quiet_NaN();
        } else {

            if (y == 0.) {
                return std::numeric_limits<double>::quiet_NaN();
            } else if (isinf_(y)) {
                return x;
            } else {
                // both are non-zero finite
                double const absx = abs_(x);
                double const absy = abs_(y); // sign of y has no effect
                double dividend = absx;
                double divisor = absy;
                bool last_bit = false;

                if (absx >= absy) {
                    int expx;
                    int expy;
                    (void)frexp_(absx, &expx);
                    (void)frexp_(absy, &expy);
                    // expx >= expy
                    divisor = ldexp(divisor, expx - expy);

                    for (; divisor > absy;) {
                        // dividend < 2*divisor
                        if (dividend >= divisor) {
                            dividend -= divisor;
                        }
                        divisor *= 0.5;
                    }

                    // dividend -> [0, 2*divisor)
                    // diviser == absy

                    last_bit = false;
                    if (dividend >= divisor) {
                        // dividend -> [divisor, 2divisor)
                        dividend -= divisor;
                        // dividend -> [0, divisor)
                        last_bit = true;
                    } else {
                        // dividend -> [0, divisor)
                        // last_bit = false;
                    }
                    // dividend -> [0, divisor)
                } else {
                    // dividend -> [0, divisor)
                }
                // diviser == absy

                double dividend_ = dividend;
                double divisor_ = divisor;

                if (0.5 * divisor * 2. != divisor) {
                    dividend_ *= 2;
                    divisor_ *= 2.;
                }

                if (dividend_ > 0.5 * divisor_) {
                    // dividend -> £¨divisor/2, divisor)
                    dividend_ -= divisor_;
                    // dividend -> £¨-divisor/2, 0)
                } else if (dividend_ < 0.5 * divisor_) {
                    // dividend -> [0, divisor/2)
                } else if (dividend_ == 0.5 * divisor_) {
                    if (last_bit) {
                        dividend_ -= divisor_;
                        // dividend -> -divisor/2
                    } else {
                        // dividend -> divisor/2
                    }
                }

                if (divisor_ != divisor) {
                    dividend_ /= 2.;
                }

                if (x < 0.) {
                    return -dividend_;
                } else {
                    return dividend_;
                }
            }
        }
    }

    inline constexpr double remainder(double x, double y) noexcept
    {
        if (std::is_constant_evaluated()) {
            return remainder_(x, y);
        } else {
            return std::fmod(x, y);
        }
    }

}


namespace qd {
    /* Computes the nearest integer to d. */
    inline constexpr double nint(double d)
    {
        // not the same as floor(d + 0.5)
        // floor(double(0.5 + 0.49999999999999994)) -> 1
        // note: 0.49999999999999994 == nextafter(0.5, 0)
        return fb::round(d);
    }

    /* Computes the truncated integer. */
    inline constexpr double aint(double d)
    {
        return (d >= 0.0) ? fb::floor(d) : fb::ceil(d);
    }
}

/* These are provided to give consistent
   interface for double with double-double and quad-double. */
inline void sincosh(double t, double& sinh_t, double& cosh_t)
{
    sinh_t = std::sinh(t);
    cosh_t = std::cosh(t);
}

inline constexpr double sqr(double t)
{
    return t * t;
}
inline constexpr double to_double(double a) { return a; }
inline constexpr int    to_int(double a) { return static_cast<int>(a); }

#endif
/*  util.h  */
#ifndef QD_UTILS_H
#define QD_UTILS_H
#include <string>

#include <iostream>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>

namespace qd {

    inline void error(const char* msg);
    /* This routine is called whenever a fatal error occurs. */
    inline void error(const char* msg)
    {
        if (msg) { std::cerr << "ERROR " << msg << std::endl; }
    }

}
void append_expn(std::string &str, int expn);


// report error by return -1
// return 0 if everything is OK
constexpr inline int read_int(char const* s, int& e)
{
    int v = 0;
    int sign = 0;
    if (*s == '-') {
        sign = -1;
        s += 1;
    } else if (*s == '+') {
        sign = '+';
        s += 1;
    }

    for (; *s; ++s) {
        char ch = *s;
        if (ch >= '0' && ch <= '9') {
            if (v < INT_MIN / 10) {
                // overflow
                return -1;
            }
            v *= 10;
            int d = *s - '0';
            v += d;
            if (v < 0) {
                // overflow
                break;
            }
        }
    }

    if (sign == -1) {
        if (v == INT_MIN) {
            e = INT_MIN;
        } else if (v < 0) {
            // true overflow
            return -1;
        } else {
            e = -v;
        }
    } else {
        e = v;
    }
    return 0;

}

// report error by return INT_MIN
constexpr inline int to_int(char const* s)
{
    int a;
    if (read_int(s, a) < 0) {
        return INT_MIN;
    }
    return a;
}
/*  util.inl.h  */
#include <cstdlib>

inline void append_expn(std::string &str, int expn) {
  int k;

  str += (expn < 0 ? '-' : '+');
  expn = std::abs(expn);

  if (expn >= 100) {
    k = (expn / 100);
    str += '0' + k;
    expn -= 100*k;
  }

  k = (expn / 10);
  str += '0' + k;
  expn -= 10*k;

  str += '0' + expn;
}

#endif


QD_CONSTEXPR bool dd_real::isnan() const { return QD_ISNAN(x[0]) || QD_ISNAN(x[1]); }
QD_CONSTEXPR bool dd_real::isfinite() const { return QD_ISFINITE(x[0]); }
QD_CONSTEXPR bool dd_real::isinf() const { return QD_ISINF(x[0]); }

#define QD_OP_TAIL(x) if(!std::is_constant_evaluated()) assert(!(x).isnan())

/* Cast to double. */
inline QD_CONSTEXPR double to_double(const dd_real& a)
{
    return a.x[0];
}

/*********** Additions ************/
/* double-double = double + double */
inline QD_CONSTEXPR dd_real dd_real::add(double a, double b) {
  double s, e;
  s = qd::two_sum(a, b, e);
  return dd_real(s, e);
}

/* double-double + double */
inline QD_CONSTEXPR dd_real operator+(const dd_real &a, double b) {
  double s1, s2;
  s1 = qd::two_sum(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/* double-double + double-double */
inline QD_CONSTEXPR dd_real dd_real::ieee_add(const dd_real &a, const dd_real &b) {
  /* This one satisfies IEEE style error bound, 
     due to K. Briggs and W. Kahan.                   */
  double s1, s2, t1, t2;

  s1 = qd::two_sum(a.x[0], b.x[0], s2);
  t1 = qd::two_sum(a.x[1], b.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

inline QD_CONSTEXPR dd_real dd_real::sloppy_add(const dd_real &a, const dd_real &b) {
  /* This is the less accurate version ... obeys Cray-style
     error bound. */
  double s, e;

  s = qd::two_sum(a.x[0], b.x[0], e);
  e += (a.x[1] + b.x[1]);
  s = qd::quick_two_sum(s, e, e);
  return dd_real(s, e);
}

inline QD_CONSTEXPR dd_real operator+(const dd_real &a, const dd_real &b) {
#ifndef QD_IEEE_ADD
  return dd_real::sloppy_add(a, b);
#else
  return dd_real::ieee_add(a, b);
#endif
}

/* double + double-double */
inline QD_CONSTEXPR dd_real operator+(double a, const dd_real &b) {
  return (b + a);
}


/*********** Self-Additions ************/
/* double-double += double */
inline QD_CONSTEXPR dd_real &dd_real::operator+=(double a) {
  double s1, s2;
  s1 = qd::two_sum(x[0], a, s2);
  s2 += x[1];
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

/* double-double += double-double */
inline QD_CONSTEXPR dd_real &dd_real::operator+=(const dd_real &a) {
#ifndef QD_IEEE_ADD
  double s, e;
  s = qd::two_sum(x[0], a.x[0], e);
  e += x[1];
  e += a.x[1];
  x[0] = qd::quick_two_sum(s, e, x[1]);
  return *this;
#else
  double s1, s2, t1, t2;
  s1 = qd::two_sum(x[0], a.x[0], s2);
  t1 = qd::two_sum(x[1], a.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
#endif
}

/*********** Subtractions ************/
/* double-double = double - double */
inline QD_CONSTEXPR dd_real dd_real::sub(double a, double b) {
  double s, e;
  s = qd::two_diff(a, b, e);
  return dd_real(s, e);
}

/* double-double - double */
inline QD_CONSTEXPR dd_real operator-(const dd_real &a, double b) {
  double s1, s2;
  s1 = qd::two_diff(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/* double-double - double-double */
inline QD_CONSTEXPR dd_real operator-(const dd_real &a, const dd_real &b) {
#ifndef QD_IEEE_ADD
  double s, e;
  s = qd::two_diff(a.x[0], b.x[0], e);
  e += a.x[1];
  e -= b.x[1];
  s = qd::quick_two_sum(s, e, e);
  return dd_real(s, e);
#else
  double s1, s2, t1, t2;
  s1 = qd::two_diff(a.x[0], b.x[0], s2);
  t1 = qd::two_diff(a.x[1], b.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
#endif
}

/* double - double-double */
inline QD_CONSTEXPR dd_real operator-(double a, const dd_real &b) {
  double s1, s2;
  s1 = qd::two_diff(a, b.x[0], s2);
  s2 -= b.x[1];
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/*********** Self-Subtractions ************/
/* double-double -= double */
inline QD_CONSTEXPR dd_real &dd_real::operator-=(double a) {
  double s1, s2;
  s1 = qd::two_diff(x[0], a, s2);
  s2 += x[1];
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

/* double-double -= double-double */
inline QD_CONSTEXPR dd_real &dd_real::operator-=(const dd_real &a) {
#ifndef QD_IEEE_ADD
  double s, e;
  s = qd::two_diff(x[0], a.x[0], e);
  e += x[1];
  e -= a.x[1];
  x[0] = qd::quick_two_sum(s, e, x[1]);
  return *this;
#else
  double s1, s2, t1, t2;
  s1 = qd::two_diff(x[0], a.x[0], s2);
  t1 = qd::two_diff(x[1], a.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
#endif
}

inline QD_CONSTEXPR dd_real dd_real::operator+() const
{
    return dd_real(+x[0], +x[1]);
}

/*********** Unary Minus ***********/
inline QD_CONSTEXPR dd_real dd_real::operator-() const {
  return dd_real(-x[0], -x[1]);
}

/*********** Multiplications ************/
/* double-double = double * double */
inline QD_CONSTEXPR dd_real dd_real::mul(double a, double b) {
  double p, e;
  p = qd::two_prod(a, b, e);
  return dd_real(p, e);
}

/* double-double * double,  where double is a power of 2. */
inline QD_CONSTEXPR dd_real mul_pwr2(const dd_real &a, double b) {
  return dd_real(a.x[0] * b, a.x[1] * b);
}

/* double-double * double */
inline QD_CONSTEXPR dd_real operator*(const dd_real &a, double b) {
  double p1, p2;

  p1 = qd::two_prod(a.x[0], b, p2);
  p2 += (a.x[1] * b);
  p1 = qd::quick_two_sum(p1, p2, p2);
  return dd_real(p1, p2);
}

/* double-double * double-double */
inline QD_CONSTEXPR dd_real operator*(const dd_real &a, const dd_real &b) {
  double p1, p2;

  p1 = qd::two_prod(a.x[0], b.x[0], p2);
  p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
  p1 = qd::quick_two_sum(p1, p2, p2);
  QD_OP_TAIL(dd_real(p1, p2));
  return dd_real(p1, p2);
}

/* double * double-double */
inline QD_CONSTEXPR dd_real operator*(double a, const dd_real &b) {
    QD_OP_TAIL(b * a);
    return (b * a);
}

/*********** Self-Multiplications ************/
/* double-double *= double */
inline QD_CONSTEXPR dd_real &dd_real::operator*=(double a) {
  double p1, p2;
  p1 = qd::two_prod(x[0], a, p2);
  p2 += x[1] * a;
  x[0] = qd::quick_two_sum(p1, p2, x[1]);
  QD_OP_TAIL(*this);
  return *this;
}

/* double-double *= double-double */
inline QD_CONSTEXPR dd_real &dd_real::operator*=(const dd_real &a) {
  double p1, p2;
  p1 = qd::two_prod(x[0], a.x[0], p2);
  p2 += a.x[1] * x[0];
  p2 += a.x[0] * x[1];
  x[0] = qd::quick_two_sum(p1, p2, x[1]);
  QD_OP_TAIL(*this);
  return *this;
}

/*********** Divisions ************/
inline QD_CONSTEXPR dd_real dd_real::div(double a, double b) {
  double q1, q2;
  double p1, p2;
  double s, e;

  q1 = a / b;

  /* Compute  a - q1 * b */
  p1 = qd::two_prod(q1, b, p2);
  s = qd::two_diff(a, p1, e);
  e -= p2;

  /* get next approximation */
  q2 = (s + e) / b;

  s = qd::quick_two_sum(q1, q2, e);

  return dd_real(s, e);
}

/* double-double / double */
inline QD_CONSTEXPR dd_real operator/(const dd_real &a, double b) {

  double q1, q2;
  double p1, p2;
  double s, e;
  dd_real r;
  
  q1 = a.x[0] / b;   /* approximate quotient. */

  /* Compute  this - q1 * d */
  p1 = qd::two_prod(q1, b, p2);
  s = qd::two_diff(a.x[0], p1, e);
  e += a.x[1];
  e -= p2;
  
  /* get next approximation. */
  q2 = (s + e) / b;

  /* renormalize */
  r.x[0] = qd::quick_two_sum(q1, q2, r.x[1]);

  return r;
}

inline QD_CONSTEXPR dd_real dd_real::sloppy_div(const dd_real &a, const dd_real &b) {
  double s1, s2;
  double q1, q2;
  dd_real r;

  q1 = a.x[0] / b.x[0];  /* approximate quotient */

  /* compute  this - q1 * dd */
  r = b * q1;
  s1 = qd::two_diff(a.x[0], r.x[0], s2);
  s2 -= r.x[1];
  s2 += a.x[1];

  /* get next approximation */
  q2 = (s1 + s2) / b.x[0];

  /* renormalize */
  r.x[0] = qd::quick_two_sum(q1, q2, r.x[1]);
  return r;
}

inline QD_CONSTEXPR dd_real dd_real::accurate_div(const dd_real &a, const dd_real &b) {
  double q1, q2, q3;
  dd_real r;

  q1 = a.x[0] / b.x[0];  /* approximate quotient */

  r = a - q1 * b;
  
  q2 = r.x[0] / b.x[0];
  r -= (q2 * b);

  q3 = r.x[0] / b.x[0];

  q1 = qd::quick_two_sum(q1, q2, q2);
  r = dd_real(q1, q2) + q3;
  return r;
}

/* double-double / double-double */
inline QD_CONSTEXPR  dd_real operator/(const dd_real &a, const dd_real &b) {
#ifdef QD_SLOPPY_DIV
  return dd_real::sloppy_div(a, b);
#else
  return dd_real::accurate_div(a, b);
#endif
}

/* double / double-double */
inline QD_CONSTEXPR dd_real operator/(double a, const dd_real &b) {
  return dd_real(a) / b;
}

inline QD_CONSTEXPR dd_real inv(const dd_real &a) {
  return 1.0 / a;
}

/*********** Self-Divisions ************/
/* double-double /= double */
inline QD_CONSTEXPR dd_real &dd_real::operator/=(double a) {
  *this = *this / a;
  return *this;
}

/* double-double /= double-double */
inline QD_CONSTEXPR dd_real &dd_real::operator/=(const dd_real &a) {
  *this = *this / a;
  return *this;
}

/*********** Squaring **********/
inline QD_CONSTEXPR dd_real sqr(const dd_real& a)
{
    double p1, p2;
    double s1, s2;
    p1 = qd::two_sqr(a.x[0], p2);
    p2 += 2.0 * a.x[0] * a.x[1];
    p2 += a.x[1] * a.x[1];
    s1 = qd::quick_two_sum(p1, p2, s2);
    return dd_real(s1, s2);
}

inline QD_CONSTEXPR dd_real dd_real::sqr(double a)
{
    double p1, p2;
    p1 = qd::two_sqr(a, p2);
    return dd_real(p1, p2);
}

/* Computes the n-th power of a double-double number.
   NOTE:  0^0 causes an error.                         */
inline QD_CONSTEXPR dd_real npwr(const dd_real& a, int n)
{

    if (n == 0) {
        if (a.is_zero()) {
            return dd_real::_nan;
        }
        return 1.0;
    }

    dd_real r = a;
    dd_real s = 1.0;
    int N = fb::abs(n);

    if (N > 1) {
        /* Use binary exponentiation */
        while (N > 0) {
            if (N % 2 == 1) {
                s *= r;
            }
            N /= 2;
            if (N > 0)
                r = sqr(r);
        }
    } else {
        s = r;
    }

    /* Compute the reciprocal if n is negative. */
    if (n < 0)
        return (1.0 / s);

    return s;
}


/********** Exponentiation **********/
inline constexpr dd_real dd_real::operator^(int n) {
  return npwr(*this, n);
}


/*********** Assignments ************/
/* double-double = double */
inline QD_CONSTEXPR dd_real &dd_real::operator=(double a) {
  x[0] = a;
  x[1] = 0.0;
  return *this;
}

/*********** Equality Comparisons ************/
/* double-double == double */
inline QD_CONSTEXPR bool operator==(const dd_real &a, double b) {
  return (a.x[0] == b && a.x[1] == 0.0);
}

/* double-double == double-double */
inline QD_CONSTEXPR bool operator==(const dd_real &a, const dd_real &b) {
  return (a.x[0] == b.x[0] && a.x[1] == b.x[1]);
}

/* double == double-double */
inline QD_CONSTEXPR bool operator==(double a, const dd_real &b) {
  return (a == b.x[0] && b.x[1] == 0.0);
}

/*********** Greater-Than Comparisons ************/
/* double-double > double */
inline QD_CONSTEXPR bool operator>(const dd_real &a, double b) {
  return (a.x[0] > b || (a.x[0] == b && a.x[1] > 0.0));
}

/* double-double > double-double */
inline QD_CONSTEXPR bool operator>(const dd_real &a, const dd_real &b) {
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] > b.x[1]));
}

/* double > double-double */
inline QD_CONSTEXPR bool operator>(double a, const dd_real &b) {
  return (a > b.x[0] || (a == b.x[0] && b.x[1] < 0.0));
}

/*********** Less-Than Comparisons ************/
/* double-double < double */
inline QD_CONSTEXPR bool operator<(const dd_real &a, double b) {
  return (a.x[0] < b || (a.x[0] == b && a.x[1] < 0.0));
}

/* double-double < double-double */
inline QD_CONSTEXPR bool operator<(const dd_real &a, const dd_real &b) {
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] < b.x[1]));
}

/* double < double-double */
inline QD_CONSTEXPR bool operator<(double a, const dd_real &b) {
  return (a < b.x[0] || (a == b.x[0] && b.x[1] > 0.0));
}

/*********** Greater-Than-Or-Equal-To Comparisons ************/
/* double-double >= double */
inline QD_CONSTEXPR bool operator>=(const dd_real &a, double b) {
  return (a.x[0] > b || (a.x[0] == b && a.x[1] >= 0.0));
}

/* double-double >= double-double */
inline QD_CONSTEXPR bool operator>=(const dd_real &a, const dd_real &b) {
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] >= b.x[1]));
}

/* double >= double-double */
inline QD_CONSTEXPR bool operator>=(double a, const dd_real &b) {
  return (b <= a);
}

/*********** Less-Than-Or-Equal-To Comparisons ************/
/* double-double <= double */
inline QD_CONSTEXPR bool operator<=(const dd_real &a, double b) {
  return (a.x[0] < b || (a.x[0] == b && a.x[1] <= 0.0));
}

/* double-double <= double-double */
inline QD_CONSTEXPR bool operator<=(const dd_real &a, const dd_real &b) {
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] <= b.x[1]));
}

/* double <= double-double */
inline QD_CONSTEXPR bool operator<=(double a, const dd_real &b) {
  return (b >= a);
}

/*********** Not-Equal-To Comparisons ************/
/* double-double != double */
inline constexpr bool operator!=(const dd_real &a, double b) {
  return (a.x[0] != b || a.x[1] != 0.0);
}

/* double-double != double-double */
inline constexpr bool operator!=(const dd_real &a, const dd_real &b) {
  return (a.x[0] != b.x[0] || a.x[1] != b.x[1]);
}

/* double != double-double */
inline constexpr bool operator!=(double a, const dd_real &b) {
  return (a != b.x[0] || b.x[1] != 0.0);
}

/*********** Micellaneous ************/
/*  this == 0 */
inline constexpr bool dd_real::is_zero() const {
  return (x[0] == 0.0);
}

/*  this == 1 */
inline constexpr bool dd_real::is_one() const {
  return (x[0] == 1.0 && x[1] == 0.0);
}

/*  this > 0 */
inline constexpr bool dd_real::is_positive() const {
  return (x[0] > 0.0);
}

/* this < 0 */
inline constexpr bool dd_real::is_negative() const {
  return (x[0] < 0.0);
}

/* Absolute value */
inline constexpr dd_real abs(const dd_real &a) {
  return (a.x[0] < 0.0) ? -a : a;
}

inline constexpr dd_real fabs(const dd_real &a) {
  return abs(a);
}

/* Cast to int. */
inline constexpr int to_int(const dd_real& a)
{
    return static_cast<int>(a.x[0]);
}

/* Reads in a double-double number from the string s. */
inline constexpr int dd_real::read(const char* s, dd_real& a)
{
    const char* p = s;
    char ch;
    int sign = 0;
    int point = -1;
    int nd = 0;
    int e = 0;
    bool done = false;
    dd_real r = 0.0;

    /* Skip any leading spaces */
    while (*p == ' ')
        p++;

    if (*p == '-') {
        sign = -1;
        ++p;
    } else if (*p == '+') {
        sign = 1;
        ++p;
    }

    while (!done && (ch = *p) != '\0') {
        if (ch >= '0' && ch <= '9') {
            int d = ch - '0';
            r *= 10.0;
            r += static_cast<double>(d);
            nd++;
        } else {

            switch (ch) {

            case '.':
                if (point >= 0)
                    return -1;
                point = nd;
                break;

            case 'E':
            case 'e':
                ++p;
                if (read_int(p, e) < 0) {
                    return -1;
                }
                done = true;
                break;

            default:
                return -1;
            }
        }

        p++;
    }

    if (point >= 0) {
        e -= (nd - point);
    }

    if (e != 0) {
        r *= (dd_real(10.0) ^ e);
    }

    a = (sign == -1) ? -r : r;
    return 0;

}

inline constexpr dd_real dd_real::read(const char* s)
{
    dd_real a;
    if (read(s, a) < 0)
        return dd_real::_nan;
    return a;
}

#endif /* _QD_DD_INLINE_H */
#endif /* _QD_DD_REAL_H */

/*  dd_math.h  */
/*
 * src/dd_real.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Contains implementation of non-inlined functions of double-double
 * package.  Inlined functions are found in dd_inline.h (in include directory).
 */
#ifndef _QD_DD_MATH_H
#define _QD_DD_MATH_H
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm> //HL: for std::max

/*  bits.h  */
/*
 * include/bits.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This file defines various routines to get / set bits of a IEEE floating
 * point number.  This is used by the library for debugging purposes.
 */

#ifndef _QD_BITS_H
#define _QD_BITS_H

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>

/* Returns the exponent of the double precision number.
   Returns INT_MIN is x is zero, and INT_MAX if x is INF or NaN. */
int get_double_expn(double x);

/* Prints 
     SIGN  EXPN  MANTISSA
   of the given double.  If x is NaN, INF, or Zero, this
   prints out the strings NaN, +/- INF, and 0.             */
void print_double_info(std::ostream &os, double x);

/*  bits.inl.h  */
/*
 * src/bits.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Defines various routines to get / set bits of a IEEE floating point
 * number.  This used by the library for debugging purposes.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <climits>



inline int get_double_expn(double x) {
  if (x == 0.0)
    return INT_MIN;
  if (QD_ISINF(x) || QD_ISNAN(x))
    return INT_MAX;

  double y = std::abs(x);
  int i = 0;
  if (y < 1.0) {
    while (y < 1.0) {
      y *= 2.0;
      i++;
    }
    return -i;
  } else if (y >= 2.0) {
    while (y >= 2.0) {
      y *= 0.5;
      i++;
    }
    return i;
  }
  return 0;
}

inline void print_double_info(std::ostream &os, double x) {
  std::streamsize old_prec = os.precision(19);
  std::ios_base::fmtflags old_flags  = os.flags();
  os << std::scientific;

  os << std::setw(27) << x << ' ';
  if (QD_ISNAN(x) || QD_ISINF(x) || (x == 0.0)) {
    os << "                                                           ";
  } else {

    x = std::abs(x);
    int expn = get_double_expn(x);
    double d = std::ldexp(1.0, expn);
    os << std::setw(5) << expn << " ";
    for (int i = 0; i < 53; i++) {
      if (x >= d) {
        x -= d;
        os << '1';
      } else
        os << '0';
      d *= 0.5;
    }

    if (x != 0.0) {
      // should not happen
      os << " +trailing stuff";
    }
  }

  os.precision(old_prec);
  os.flags(old_flags);
}

#endif  /* _QD_BITS_H */

/*  double_math.h  */
#ifndef QD_FLOATING_MATH_H
#define QD_FLOATING_MATH_H


namespace fb {

#if 0
    inline constexpr double fma__(double a, double b, double c) noexcept
    {
        if (isnan_(a) || isnan_(b) || isnan_(c)) {
            return c;
        } else if (isinf_(a)) {
            if (b == 0.) {
                return _d_nan;
            } else if (isinf_(c)) {
                if ((signbit(a) == signbit(b)) != signbit(c)) {
                    return _d_nan;
                } else {
                    return c;
                }
            }
        } else if (isinf_(b)) {
            if (a == 0.) {
                return _d_nan;
            } else if (isinf_(c)) {
                if ((signbit(a) == signbit(b)) != signbit(c)) {
                    return _d_nan;
                } else {
                    return c;
                }
            }
        } else if (isinf_(c)) {
            return c;
        } else {
            // a,b,c are all finite
            // calculate a*b+c

            bool absign = signbit(a) == signbit(b);
            bool sign = absign == signbit(c);
            int aexp, bexp, cexp;
            double a_ = frexp(a, &aexp);
            double b_ = frexp(b, &bexp);
            double c_ = frexp(c, &cexp);

            uint64_t af = frac_part(a_);
            uint64_t bf = frac_part(b_);
            uint64_t cf_ = frac_part(c_);
            af |= uint64_t(1) << EXP_OFFSET;
            bf |= uint64_t(1) << EXP_OFFSET;
            cf_ |= uint64_t(1) << EXP_OFFSET;
            Uint128 abf = Uint128::mul(af, bf);
            uint64_t cf = cf_;

            int abfexp = aexp + bexp - 52 - 52;
            // a*b = abf * 2^{ abfexp }
            int cfexp = cexp - 52;
            // c = cf * 2^{ cfexp }

            int coff = cfexp - abfexp;

            Uint128 ans = 0;
            int ansexp = 0;
            if (sign) {
                if (coff <= -53) {
                    ans = abf;
                    ansexp = abfexp;
                } else if (coff >= 107) { // 53bits * 53*bits = 107 bits
                    ans = cf;
                    ansexp = cfexp;
                } else {

                }

            }

        }
    }
#endif

    inline constexpr double fdim_(double x, double y)
    {
        if (x > y) {
            return x - y;
        } else if (x <= y) {
            return +0.;
        } else {
            return _d_nan;
        }
    }

    inline constexpr double nan_(char const* s)
    {
        if (*s == '\0') {
            return _d_nan;
        }
        uint64_t v = 0;
        for (; *s; ++s) {
            char ch = *s;
            if (ch >= '0' && ch <= '9') {
                v *= 10;
                int d = *s - '0';
                v += d;
            }
        }
        // TODO: negative? overflow?

        constexpr int exp = exp_part(_d_nan);
        return make_double(v, exp, 0);
    }

    inline constexpr double log1p_Taylor(double x)
    {
        double x_frac = x;
        double z1 = x_frac;
        double z2 = 2. + x_frac;
        double z3 = z1 / z2;
        double z6 = sqr_(z3);


        double t = 0;
        t = t * z6 + 1 / 13.;
        t = t * z6 + 1 / 11;
        t = t * z6 + 1 / 9.;
        t = t * z6 + 1 / 7.;
        t = t * z6 + 1 / 5.;
        t = t * z6 + 1 / 3.;
        t = t * z6 + 1 / 1.;

        double r = 2 * z3 * t;
        return r;
    }

    //0.5 <= x < 1
    inline constexpr double log1p_2(double x)
    {
        constexpr double _d_ln_1d5 = 0.4054651081081643819780131154643491365719904234624941976140143241;

        constexpr double _d_ln_table[] = {
            0.4054651081081643819780131154643491365719904234624941976140143241,
            0.2231435513142097557662950903098345033746010855480072136712878724,
             0.1177830356563834545387941094705217050684807125647331411073486387,
             0.0606246218164348425806061320404202632862024751447237708145176999,
             0.0307716586667536883710282075967721640916967399588903563498619953,
             0.0155041865359652541508540460424468358778684928671933136076133451,
             0.0077821404420549489474629000611367636781258021825180880816195321,
             0.0038986404156573230139373430958429070107237541049028050776797502,
             0.0019512201312617494396740495318415385003497255250798866559231518,
             0.0009760859730554588959608249080171866726118343337845362377585982,
        };

        double lnx_ = 0.;
        double x_frac = x;

        if (x_frac > 0) {
            if (x_frac >= 0.5) {
                // we want (1 + x_frac) / 1.5
                // (1.5 + x_frac - 0.5) / 1.5
                // 1 + (x_frac - 0.5)/1.5
                x_frac = (x_frac - 0.5) / 1.5;
                // x_frac < 1 -> x_frac < 0.5
                lnx_ += _d_ln_1d5;
            } else {
                // x_frac < 0.5
            }
            // x_frac < 0.5

            bool _d_ln_table_v[std::size(_d_ln_table)] = {};
            double frac_test = 1 / 4.;
            for (size_t i = 1; i < std::size(_d_ln_table); ++i) {
                if (x_frac >= frac_test) {
                    x_frac = (x_frac - frac_test) / (1 + frac_test);
                    _d_ln_table_v[i] = true;
                }
                frac_test /= 2;
            }

            double t2 = 0;
            for (size_t i = std::size(_d_ln_table) - 1; i >= 1; --i) {
                t2 += _d_ln_table_v[i] ? _d_ln_table[i] : 0;
            }
            lnx_ += t2;

        } else {

            bool _d_ln_table_v[std::size(_d_ln_table)] = {};

            double frac_test = 0.5;
            for (size_t i = 0; i < std::size(_d_ln_table); ++i) {

                // 0 > x >= -test;
                // (1+x)*(1+frac_test) - 1
                // x + (1+ frac_test) * x
                double x2 = frac_test + (1 + frac_test) * x_frac;
                // x2 = frac_test + x + frac_test*x
                // x2 >= frac_test - frac_test - frac_test*frac_test
                // x2 >= - frac_test*frac_test >= -0.5 frac_test
                if (x2 <= 0) {
                    _d_ln_table_v[i] = true;
                    x_frac = x2;
                }
                frac_test /= 2;
            }

            double t2 = 0;
            for (size_t i = std::size(_d_ln_table); i >= 1; --i) {
                t2 += _d_ln_table_v[i - 1] ? _d_ln_table[i - 1] : 0;
            }
            lnx_ -= t2;
        }

        double fix = log1p_Taylor(x_frac);
        lnx_ += fix;

        return lnx_;
    }


    inline constexpr double log_(double x)
    {
        if (isnan(x)) {
            return _d_nan;
        } else if (x == 0.) {
            return copysign_(_infinity, x);
        } else if (x == 1) {
            return 0.;
        } else if (x < 0.) {
            return _d_nan;
        } else if (isinf(x)) {
            return x;
        } else {
            //positive finite

            if (x >= 0.5 && x < 2.) {
                return log1p_2(x - 1);
            } else {
                int exp_;
                double x_ = frexp_(x, &exp_);
                x_ *= 2.;
                exp_ -= 1;

                // x = x_ 2^exp_ 
                // 1. <= x_ < 2
                // log(x) = exp_ * log(x_)
                // log(x) = log(x_) + exp_*log(2)

                double r = log1p_2(x_ - 1.);
                r += exp_ * _d_ln2;
                return r;

            }
        }
    }

    inline constexpr double log(double x) noexcept
    {
        if (std::is_constant_evaluated()) {
            return log_(x);
        } else {
            return std::log(x);
        }
    }

    inline constexpr double log1p_(double x)
    {
        if (isnan_(x)) {
            return _d_nan;
        } else if (x == 0.) {
            return x;
        } else if (x < -1.) {
            return _d_nan;
        } else if (isinf_(x)) {
            return _infinity;
        } else {
            //positive finite

            if (x < -0.5) {
                return log_(1. + x);
            } else if (x > 1.) {
                return log_(1. + x);
            } else {
                return log1p_2(x);
            }
        }
    }

    inline constexpr double log1p(double x)
    {
        if (std::is_constant_evaluated()) {
            return log1p_(x);
        } else {
            return std::log1p(x);
        }
    }

    inline constexpr double log10_(double x) noexcept
    {
        return log_(x) / _d_ln10;
    }

    inline constexpr double log2_(double x) noexcept
    {
        int a;
        double x_ = frexp_(x, &a);
        return log(x_) / _d_ln2 + a;
    }
#if 0
    inline bool self_add_car(uint64_t& a, uint64_t b)
    {
        a += b;
        return a < b;
    }
    bool get_double(uint64_t bits, uint64_t remain, double& r)
    {
        int lb = leading_bit(bits);
        if (lb >= 0) {
            uint64_t ur = bits << (52 - lb);
            ur |= (remain >> 12) + lb;
            r = ldexp(double(ur), -53 + lb);
            return true;
        } else {
            return false;
        }
    }
#endif

    // a mod pi/2
    inline constexpr dd_real mod_pio2_dd(double a, int& n)
    {
        if (abs_(a) < _d_pi / 4) {
            n = 0;
            return a;
        }

        bool sign = false;
        if (a < 0) {
            a = -a;
            sign = true;
        }
        int m = 0;



        // 1152 = 18 * 64 bits for 2/pi
        bool const _2opi[] = {
            1,0,1,0,0,0,1,0,1,1,1,1,1,0,0,1,
1,0,0,0,0,0,1,1,0,1,1,0,1,1,1,0,
0,1,0,0,1,1,1,0,0,1,0,0,0,1,0,0,
0,0,0,1,0,1,0,1,0,0,1,0,1,0,0,1,
1,1,1,1,1,1,0,0,0,0,1,0,0,1,1,1,
0,1,0,1,0,1,1,1,1,1,0,1,0,0,0,1,
1,1,1,1,0,1,0,1,0,0,1,1,0,1,0,0,
1,1,0,1,1,1,0,1,1,1,0,0,0,0,0,0,
1,1,0,1,1,0,1,1,0,1,1,0,0,0,1,0,
1,0,0,1,0,1,0,1,1,0,0,1,1,0,0,1,
0,0,1,1,1,1,0,0,0,1,0,0,0,0,1,1,
1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,
1,1,1,1,1,1,1,0,0,1,0,1,0,0,0,1,
0,1,1,0,0,0,1,1,1,0,1,0,1,0,1,1,
1,1,0,1,1,1,1,0,1,0,1,1,1,0,1,1,
1,1,0,0,0,1,0,1,0,1,1,0,0,0,0,1,
1,0,1,1,0,1,1,1,0,0,1,0,0,1,0,0,
0,1,1,0,1,1,1,0,0,0,1,1,1,0,1,0,
0,1,0,0,0,0,1,0,0,1,0,0,1,1,0,1,
1,1,0,1,0,0,1,0,1,1,1,0,0,0,0,0,
0,0,0,0,0,1,1,0,0,1,0,0,1,0,0,1,
0,0,1,0,1,1,1,0,1,1,1,0,1,0,1,0,
0,0,0,0,1,0,0,1,1,1,0,1,0,0,0,1,
1,0,0,1,0,0,1,0,0,0,0,1,1,1,0,0,
1,1,1,1,1,1,1,0,0,0,0,1,1,1,0,1,
1,1,1,0,1,0,1,1,0,0,0,1,1,1,0,0,
1,0,1,1,0,0,0,1,0,0,1,0,1,0,0,1,
1,0,1,0,0,1,1,1,0,0,1,1,1,1,1,0,
1,1,1,0,1,0,0,0,1,0,0,0,0,0,1,0,
0,0,1,1,0,1,0,1,1,1,1,1,0,1,0,1,
0,0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,
0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,0,
1,1,1,0,1,0,0,1,1,0,0,1,1,1,0,0,
0,1,1,1,0,0,0,0,0,0,1,0,0,1,1,0,
1,0,1,1,0,1,0,0,0,1,0,1,1,1,1,1,
0,1,1,1,1,1,1,0,0,1,0,0,0,0,0,1,
0,0,1,1,1,0,0,1,1,0,0,1,0,0,0,1,
1,1,0,1,0,1,1,0,0,0,1,1,1,0,0,1,
1,0,0,0,0,0,1,1,0,1,0,1,0,0,1,1,
0,0,1,1,1,0,0,1,1,1,1,1,0,1,0,0,
1,0,0,1,1,1,0,0,1,0,0,0,0,1,0,0,
0,1,0,1,1,1,1,1,1,0,0,0,1,0,1,1,
1,0,1,1,1,1,0,1,1,1,1,1,1,0,0,1,
0,0,1,0,1,0,0,0,0,0,1,1,1,0,1,1,
0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,
1,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,
1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,1,
1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,
1,1,1,0,1,1,1,1,0,0,1,0,1,1,1,1,
0,0,0,1,0,0,0,1,1,0,0,0,1,0,1,1,
0,1,0,1,1,0,1,0,0,0,0,0,1,0,1,0,
0,1,1,0,1,1,0,1,0,0,0,1,1,1,1,1,
0,1,1,0,1,1,0,1,0,0,1,1,0,1,1,0,
0,1,1,1,1,1,1,0,1,1,0,0,1,1,1,1,
0,0,1,0,0,1,1,1,1,1,0,0,1,0,1,1,
0,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,
0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,0,
0,0,1,1,1,1,1,1,0,1,1,0,0,1,1,0,
1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,
1,1,1,0,1,0,1,0,0,0,1,0,1,1,0,1,
0,1,1,1,0,1,0,1,0,0,1,0,0,1,1,1,
1,0,1,1,1,0,1,0,1,1,0,0,0,1,1,1,
1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,
1,1,1,1,0,0,0,1,0,1,1,1,1,0,1,1,
0,0,1,1,1,1,0,1,0,0,0,0,0,1,1,1,
0,0,1,1,1,0,0,1,1,1,1,1,0,1,1,1,
1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,0,
1,0,0,1,0,0,1,0,1,1,1,0,1,0,1,0,
0,1,1,0,1,0,1,1,1,1,1,1,1,0,1,1,
0,1,0,1,1,1,1,1,1,0,1,1,0,0,0,1,
0,0,0,1,1,1,1,1,1,0,0,0,1,1,0,1,
0,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,
0,1,0,1,0,1,1,0,0,0,0,0,0,0,1,1,
0,0,1,1,0,0,0,0,0,1,0,0,0,1,1,0,
1,1,1,1,1,1,0,0,0,1,1,1,1,0,1,1,
        };

        // calculate a * _2opi to [2 bits . 121 bits] is enough
        // of course, for abs(a) >= pi/4, if a <= pi/4 deal with it specailly
        // see ARGUMENT REDUCTION FOR HUGE ARGUMENTS: Good to the Last Bit

        constexpr int B = 20;
        uint64_t u = frac_part(a) | (uint64_t(1) << 52);
        int exp_ = exp_part(a);
        // a = u * 2^(exp_ - 52)
        // u = 2^[0,52]
        // _2opi = 2^{-i-1}
        // a * _2opi = 2^{-i-1} 2^[0,52] 2^(exp_ - 52)
        // = 2^  { exp_ - i - 53 + [0,52] }
        // exp_ - i - 53 + [0,52] = [1, -121-B]
        // i =  exp_ - [1, -121-B] + [0, 52] - 53
        // i =  exp_ + [-54, 120+B]
        // 175+B bits

        int buffer[175 + 53 + B] = {};

        for (int i_ = 0; i_ < 175 + B; ++i_) {
            int i = exp_ + 120 + B - i_;
            bool bit = 0;
            if (i >= 0) {
                bit = _2opi[i];
            }
            if (bit) {
                int car = 0;
                for (int j = 0; j < 53; ++j) {
                    buffer[j + i_] += ((u >> j) & 1) + car;
                    if (buffer[j + i_] > 1) {
                        buffer[j + i_] -= 2;
                        car = 1;
                    } else {
                        car = 0;
                    }
                }
                for (int j = 53; (j + i_) < (int)std::size(buffer) && car; ++j) {
                    buffer[j + i_] += car;
                    if (buffer[j + i_] > 1) {
                        buffer[j + i_] -= 2;
                        car = 1;
                    } else {
                        car = 0;
                    }
                }
            }
        }


        // buffer factor = 2^{exp_ - 52 + -(exp_ + 120 + B + 1)} = 2^{-173-B}
        int car = 0;
        if (buffer[172 + B]) { // frac >=0.5

            // frac' = 1 - frac
            int car_ = 0;
            for (int i = 0; i <= 172 + B; ++i) {
                int v = car_ - buffer[i];

                if (v < 0) {
                    v += 2;
                    car_ = -1;
                } else {
                    car_ = 0;
                }
                buffer[i] = v;
            }
            car = -1;
        }

        int int_part = buffer[173 + B] + buffer[174 + B] * 2;

        // (n+frac) = (n+(1-frac')) = n+1 - frac'
        int_part += -car;
        if (int_part >= 4) int_part -= 4;
        m = int_part;

        // fraction
        dd_real r_ = 0.;
        int expn = 0;
        for (int i = 172 + B; i >= 0; --i) {
            expn -= 1;
            if (buffer[i]) {
                uint64_t u1 = 0, u2 = 0;
                for (int j = 0; j < 53; ++j, --i) {
                    u1 <<= 1;
                    u1 |= (i >= 0) ? buffer[i] : 0;
                }
                for (int j = 0; j < 53; ++j, --i) {
                    u2 <<= 1;
                    u2 |= (i >= 0) ? buffer[i] : 0;
                }
                double r1 = ldexp_(double(u1), expn - 52);
                double r2 = ldexp_(double(u2), expn - 53 - 52);
                r_ = dd_real::add(r1, r2);
                break;
            }
        }
        r_ = dd_real::_pi / 2. * r_;
        if (car) r_ = -r_;

        // fix for sign
        if (sign) {
            r_ = -r_;
            if (m != 0) {
                m = 4 - m;
            }
        }

        n = m;
        return r_;
    }

    // a mod pi/2
    inline constexpr double mod_pio2(double a, int& n)
    {
        return to_double(mod_pio2_dd(a, n));
    }

    inline constexpr dd_real sin_table[] = {
    dd_real::read("0"),
    dd_real::read("0.124674733385227689957442708712108"),
    dd_real::read("0.247403959254522929596848704849389"),
    dd_real::read("0.366272529086047561372909351716264"),
    dd_real::read("0.479425538604203000273287935215571"),
    dd_real::read("0.585097272940462154805399314150080"),
    dd_real::read("0.681638760023334166733241952779893"),
    dd_real::read("0.767543502236027039634575467054539"),
    };

    inline constexpr dd_real cos_table[] = {
    dd_real::read("1"),
    dd_real::read("0.992197667229329053149096907788250"),
    dd_real::read("0.968912421710644784144595449494189"),
    dd_real::read("0.930507621912314291149476792229555"),
    dd_real::read("0.877582561890372716116281582603829"),
    dd_real::read("0.810963119505217902189534803941080"),
    dd_real::read("0.731688868873820886311838753000084"),
    dd_real::read("0.640996858163325130356556622796034"),
    };

    inline constexpr double phi_table[] = {
        0, 1 / 8.,2 / 8.,3 / 8.,4 / 8., 5 / 8., 6 / 8.,7 / 8.,
    };

    inline constexpr size_t get_nearest_phi(double x)
    {
        size_t rnd = round_(x * 8.);
        return rnd <= 7 ? rnd : 7;
    }

    inline constexpr dd_real sin_Taylor(dd_real x)
    {
        // x - 1/3! x^3 + 1/5! x^7
        // x(1 - 1/(3 2) x^2 (1 - 1/£¨4 5)x^2(1 - 1/(6 7)x^2)))

        dd_real const t = x * x;
        dd_real v = 1.;
        v = 1. - t * v / (18. * 19);
        v = 1. - t * v / (16. * 17);
        v = 1. - t * v / (14. * 15);
        v = 1. - t * v / (12. * 13);
        v = 1. - t * v / (10. * 11);
        v = 1. - t * v / (8. * 9);
        v = 1. - t * v / (6. * 7);
        v = 1. - t * v / (4. * 5);
        v = 1. - t * v / (2. * 3);
        v = v * x;
        return v;
    }


    inline constexpr dd_real cos_Taylor(dd_real x)
    {
        // 1 - 1/2! x^2 + 1/4! x^4
        // 1 - 1/(1 2) x^2 (1 - 1/(3 4)x^2 (1 - 1/(5 6) x^2))

        dd_real const t = x * x;
        dd_real v = 1.;
        v = 1. - t * v / (19. * 20);
        v = 1. - t * v / (17. * 18);
        v = 1. - t * v / (15. * 16);
        v = 1. - t * v / (13. * 14);
        v = 1. - t * v / (11. * 12);
        v = 1. - t * v / (9. * 10);
        v = 1. - t * v / (7. * 8);
        v = 1. - t * v / (5. * 6);
        v = 1. - t * v / (3. * 4);
        v = 1. - t * v / (1. * 2);
        return v;
    }


    inline constexpr dd_real sin_can(dd_real x)
    {
        bool sign = x < 0.;
        if(x < 0.) x = -x;
        size_t idx = get_nearest_phi(to_double(x));
        x = x - phi_table[idx];
        dd_real s = sin_Taylor(x);
        dd_real c = cos_Taylor(x);
        dd_real sb = sin_table[idx];
        dd_real cb = cos_table[idx];
        dd_real v = c * sb + s * cb;
        return sign ? -v : v;
    }

    inline constexpr dd_real cos_can(dd_real x)
    {
        if (x < 0.) x = -x;

        size_t idx = get_nearest_phi(to_double(x));
        x = x - phi_table[idx];
        dd_real s = sin_Taylor(x);
        dd_real c = cos_Taylor(x);
        dd_real sb = sin_table[idx];
        dd_real cb = cos_table[idx];
        dd_real v = c * cb - s * sb;
        return v;
    }

    inline constexpr double cos_(double x)
    {
        if (x == 0. || isnan_(x)) {
            return 1.;
        } else if (isinf_(x)) {
            return _d_nan;
        } else {
            // finite
            int n;
            double a = mod_pio2(x, n);
            double v = 0.;
            switch (n) {
            case 0:
                v = to_double(cos_can(a));
                break;
            case 1:
                v = -to_double(sin_can(a));
                break;
            case 2:
                v = -to_double(cos_can(a));
                break;
            case 3:
                v = to_double(sin_can(a));
                break;
            }
            return v;
        }
    }
    inline constexpr double sin_(double x)
    {
        if (x == 0. || isnan_(x)) {
            return x;
        } else if (isinf_(x)) {
            return _d_nan;
        } else {
            // finite

            int n;
            double a = mod_pio2(x, n);
            double v = 0.;
            switch (n) {
            case 0:
                v = to_double(sin_can(a));
                break;
            case 1:
                v = to_double(cos_can(a));
                break;
            case 2:
                v = -to_double(sin_can(a));
                break;
            case 3:
                v = -to_double(cos_can(a));
                break;
            }
            return v;


        }
    }

    inline constexpr double sin(double x) noexcept
    {
        if (std::is_constant_evaluated()) {
            return sin_(x);
        } else {
            return std::sin(x);
        }
    }

    inline constexpr double cos(double x) noexcept
    {
        if (std::is_constant_evaluated()) {
            return cos_(x);
        } else {
            return std::cos(x);
        }
    }

    inline constexpr double tan_(double x) noexcept
    {
        if (x == 0. || isnan_(x)) {
            return x;
        } else if (isinf_(x)) {
            return _d_nan;
        } else {
            // finite

            int n;
            double a = mod_pio2(x, n);
            double v = 0.;
            switch (n) {
            case 0:
                v = to_double(sin_can(a) / cos_can(a));
                break;
            case 1:
                v = -to_double(cos_can(a) / sin_can(a));
                break;
            case 2:
                v = to_double(sin_can(a) / cos_can(a));
                break;
            case 3:
                v = -to_double(cos_can(a) / sin_can(a));
                break;
            }
            return v;


        }
    }

    inline constexpr double tan(double x) noexcept
    {
        if (std::is_constant_evaluated()) {
            return tan_(x);
        } else {
            return std::cos(x);
        }
    }

    // error 3E-7 a^11
    inline constexpr double expm1_Taylor(double a) noexcept
    {
        // x ( 1 + 1/2 x (1 + 1/3 x(1 + 1/4 x)))
        double v = 1;
        v = 1 + 1 / 10. * a * v;
        v = 1 + 1 / 9. * a * v;
        v = 1 + 1 / 8. * a * v;
        v = 1 + 1 / 7. * a * v;
        v = 1 + 1 / 6. * a * v;
        v = 1 + 1 / 5. * a * v;
        v = 1 + 1 / 4. * a * v;
        v = 1 + 1 / 3. * a * v;
        v = 1 + 1 / 2. * a * v;
        v *= a;
        return v;
    }

#define _FB_SPLITTER 134217729.0               // = 2^27 + 1
#define _FB_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

    inline constexpr void split_double(double a, double& hi, double& lo)
    {
        double temp;
        if (a > _FB_SPLIT_THRESH || a < -_FB_SPLIT_THRESH) {
            a *= 3.7252902984619140625e-09;  // 2^-28
            temp = _FB_SPLITTER * a;
            hi = temp - (temp - a);
            lo = a - hi;
            hi *= 268435456.0;          // 2^28
            lo *= 268435456.0;          // 2^28
        } else {
            temp = _FB_SPLITTER * a;
            hi = temp - (temp - a);
            lo = a - hi;
        }
    }

    /* Computes fl(a+b) and err(a+b).  */
    inline constexpr double two_sum(double a, double b, double& err)
    {
        double s = a + b;
        double bb = s - a;
        err = (a - (s - bb)) + (b - bb);
        return s;
    }

    /* Computes fl(a*a) and err(a*a).  */
    inline constexpr double two_sqr(double a, double& err)
    {
        double hi, lo;
        double q = a * a;
        split_double(a, hi, lo);
        err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
        return q;
    }

    /* Computes fl(a*b) and err(a*b).  */
    inline constexpr double two_prod(double a, double b, double& err)
    {
        double a_hi, a_lo, b_hi, b_lo;
        double p = a * b;
        split_double(a, a_hi, a_lo);
        split_double(b, b_hi, b_lo);
        err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
        return p;
    }

    struct value_error {
        double value;
        double error;


        constexpr value_error mul2() const
        {
            value_error ans;
            ans.value = 2 * value;
            ans.error = 2 * error;
            return ans;
        }

    };

    constexpr double sqr(double x)
    {
        return x * x;
    }

    constexpr value_error sqr(value_error x)
    {
        value_error ans;
        ans.value = two_sqr(x.value, ans.error);
        ans.error += 2 * x.value * x.error + x.error * x.error;
        two_sum(ans.value, ans.error, ans.error);
        return ans;
    }

    constexpr value_error operator+(value_error x, value_error y)
    {
        value_error ans;
        ans.value = two_sqr(x.value, ans.error);
        ans.error += 2 * x.value * x.error + x.error * x.error;
        two_sum(ans.value, ans.error, ans.error);
        return ans;
    }

    inline constexpr double expm1_n(double a, int& n, double& err) noexcept
    {
        err = 0.;
        // exp(x + n ln2) = exp(x) 2^n
        double nd = floor_(a / _d_ln2);
        n = int(nd);

        double p1, p2, p3 = _d_ln2_t1;
        split_double(_d_ln2, p1, p2);
        double x = a;
        double dx1 = 0;
        double dx2 = 0;
        x = two_sum(x, n * -p1, dx1);
        x = two_sum(x, n * -(p2 + p3), dx2);
        double dx = dx1 + dx2;
        x = two_sum(x, dx, dx);




#if 0
        // [0, ln2]
        // exp(x) = exp(x-a)exp(a) = (1 + expm1(x-a))exp(a)
        // = exp(a) + expm1(x-a)*exp(a)
        // don't forget: dx*exp(x-a)exp(a)
        // make a table, result is more accuracy
        constexpr double a_table[] = {
            0 / 16., 1 / 16., 2 / 16., 3 / 16., 4 / 16.,
        5 / 16., 6 / 16.,7 / 16.,8 / 16.,9 / 16.,10 / 16.,11 / 16. };
        constexpr double exp_table[] = {
            1.0,
            1.0644944589178594295633905946428896731007254436493533015193075106,
            1.1331484530668263168290072278117938725655031317451816259128200360,
            1.2062302494209807106555860104464335480403936461999703807388699348,
            1.2840254166877414840734205680624364583362808652814630892175072968,
            1.3668379411737963628387567727212086721727332944308111731505954490,
            1.4549914146182013360537936919875185083468420209644156811952413281,
            1.5488302986341330979985519845954923375583036629258105734128604976,
            1.6487212707001281468486507878141635716537761007101480115750793116,
            1.7550546569602985572440470365989676887382375302457485300516127462,
            1.8682459574322224065018356201881044531149722837225540862147663759,
            1.9887374695822918311174773496469253668482551764105723262843912825,
        };
        constexpr double expm1_table[] = {
            0.0,
            0.0644944589178594295633905946428896731007254436493533015193075106,
            0.1331484530668263168290072278117938725655031317451816259128200360,
            0.2062302494209807106555860104464335480403936461999703807388699348,
            0.2840254166877414840734205680624364583362808652814630892175072968,
            0.3668379411737963628387567727212086721727332944308111731505954490,
            0.4549914146182013360537936919875185083468420209644156811952413281,
            0.5488302986341330979985519845954923375583036629258105734128604976,
            0.6487212707001281468486507878141635716537761007101480115750793116,
            0.7550546569602985572440470365989676887382375302457485300516127462,
            0.8682459574322224065018356201881044531149722837225540862147663759,
            0.9887374695822918311174773496469253668482551764105723262843912825,
        };
        double nd2 = floor_(16 * x);
        int n2 = (int)nd2;
        x -= nd2 / 16.;
        double expm1 = expm1_Taylor(x);

        double err1, err2;
        double t = two_prod(expm1, exp_table[n2], err1);
        expm1 = two_sum(expm1_table[n2], t, err2);
        expm1 = two_sum(expm1, err1 + err2 + dx * (1 + expm1), err);
#else
        // exp(x + dx) = exp(x)(1+dx)= (1 + expm1)(1+dx)
        // = 1 + expm1 + dx + expm1*dx
        // exp(x/k * k) = exp(x/k)^k
        // expm1 = (1+expm1)^k - 1
        // k = 2: expm1 = 2*expm1 + expm1^2

        double expm1 = expm1_Taylor(x / 16.);
        expm1 = 2 * expm1 + sqr(expm1);
        expm1 = 2 * expm1 + sqr(expm1);
        expm1 = 2 * expm1 + sqr(expm1);
        expm1 = 2 * expm1 + sqr(expm1);
        err += (1 + expm1) * dx;
        expm1 = two_sum(expm1, err, err);
#endif

        return expm1;
    }

    inline constexpr double exp_(double a) noexcept
    {
        int n;
        double err1;
        double expm1 = expm1_n(a, n, err1);
        expm1 = expm1 + err1;
        double expx = 1 + expm1;
        expx = ldexp_(expx, n);
        return expx;
    }

    inline constexpr double exp2_(double a) noexcept
    {
        if (a == 0.) return 1;
        else if (isinf_(a) && a < 0.) {
            return 0.;
        } else if (isinf(a) && a > 0.) {
            return a;
        } else if (isnan(a)) {
            return _d_nan;
        } else {
            double frac = a - floor_(a);
            if (floor_(a) >= INT_MAX) {
                return _infinity;
            } else if (floor_(a) <= INT_MIN) {
                return 0.;
            }
            double int_ = (int)floor_(a);
            // 2^ frac = exp(frac/ln(2))
            frac /= _d_ln2;
            double base = exp_(frac);
            base = ldexp(base, int_);
            return base;
        }
    }

    inline constexpr double expm1_(double a) noexcept
    {
        int n;
        double err1;

        double expm1 = expm1_n(a, n, err1);

        if (n == 0) {
            expm1 += err1;
            return expm1;
        } else {
            double expx = two_sum(1, expm1, err1);
            expx = ldexp(expx, n);
            err1 = ldexp(err1, n);
            expx -= 1;
            expx += err1;
            return expx;
        }
    }

    inline constexpr double pow_(double base, double n)
    {
        // b^n = e^{ln(b)n}
        return exp_(log_(base) * n);
    }


    inline constexpr double cbrt_(double x)
    {
        // x = x_ 2^{3m+n} where 0 <= n < 3
        // cbrt(x) = cbrt(x_ 2^n) 2^(m)
        int exp = 0;
        double x_ = frexp_(x, &exp);
        int n = exp % 3;
        int m = exp / 3;

        for (; n; --n) {
            x_ *= 2;
        }

        // cbrt(x_)
        // [0.5, 2)
        // f'(r) = 3r^2
        // r' = r - (f(x) - x)/f'(r)
        // r' = r - (r*r*r - x)/(3 r^2)
        // r' = (x/r^2 + 2 r) / 3;
        double r = 1.;
        for (int i = 0; i < 8; ++i) {
            r = 2 / 3. * r + 1 / 3. * x_ / (r * r);
        }
#if 0
        // now we need to calculate the r*r*r - x as precision as possible
        double err, err1, err2;
        double r2 = two_prod(r, r, err);
        err *= r;
        double r3 = two_prod(r2, r, err1);
        err += err1;
        double diff = two_sum(r3, -x_, err2);
        err += err2;
        diff = two_sum(diff, err, err);
        r -= diff / (3 * r * r);
#else
        dd_real r_ = r;
        r_ = 2. / dd_real(3.) * r_ + 1. / dd_real(3.) * x_ / (r_ * r_);
        r = to_double(r_);
#endif

        r = ldexp(r, m);
        return r;
    }

    inline constexpr void swap_(double& x, double& y)
    {
        double t = x;
        x = y;
        y = t;
    }

    inline constexpr double hypot_(double x, double y)
    {
        x = abs_(x);
        y = abs_(y);

        if (isinf_(x)) {
            return x;
        } else if (isinf_(y)) {
            return y;
        } else if (x == 0.) {
            return y;
        } else if (y == 0.) {
            return x;
        } else if (isnan_(x) || isnan_(y)) {
            return _d_nan;
        }

        // x, y are all non-zero finite

        if (x < y) {
            swap_(x, y);
        }
        // x >= y

        int exp;
        x = frexp_(x, &exp);
        y = ldexp_(y, -exp);

#if 0
        double err1, err2, err3, err4, err5;
        double x2 = two_prod(x, x, err1);
        double y2 = two_prod(y, y, err2);
        double z2 = two_sum(x2, y2, err3);
        double r = sqrt_(z2);
        double r2 = two_prod(r, r, err4);
        err5 = z2 - r2;
        // sqrt(x^2 + y^2) = sqrt(x2 + err1 + y2 + err2)
        // = sqrt(z2 + err3 + err1 + err2)
        // = sqrt(r2 + err5 + err3 + err1 + err2)
        // = sqrt(r^2 - err4 + err5 + err3 + err1 + err2)
        // sqrt(x^2 + y^2) = r + 0.5 (sum err) / r
        r += 0.5 * (err1 + err2 + err3 - err4 + err5) / r;
#else

        dd_real dd_x = x;
        dd_real dd_y = y;
        dd_real r2 = dd_x * dd_x + dd_y * dd_y;
        dd_real dd_r = sqrt_(x * x + y * y);
        // Newton's method to sqrt
        dd_r = 0.5 * (r2 / dd_r + dd_r);
        dd_r = 0.5 * (r2 / dd_r + dd_r);
        double r = to_double(dd_r);

#endif
        r = ldexp_(r, exp);
        return r;
    }

    inline constexpr double hypot_(double x, double y, double z)
    {
        x = abs_(x);
        y = abs_(y);
        z = abs_(z);

        if (isinf_(x) || isinf_(y) || isinf_(z)) {
            return _infinity;
        } else if (isnan_(x) || isnan_(y) || isnan_(z)) {
            return _d_nan;
        }

        // x, y, z are all finite
        if (x < y) {
            swap_(x, y);
        }
        if (x < z) {
            swap_(x, z);
        }

        // x >= y, z
        if (x == 0.) {
            return 0.;
        }

        int exp;
        x = frexp_(x, &exp);
        y = ldexp_(y, -exp);
        z = ldexp_(z, -exp);

        dd_real dd_x = x;
        dd_real dd_y = y;
        dd_real dd_z = z;
        dd_real r2 = dd_x * dd_x + dd_y * dd_y + dd_z * dd_z;
        dd_real dd_r = sqrt_(x * x + y * y + z * z);
        // Newton's method to sqrt
        dd_r = 0.5 * (r2 / dd_r + dd_r);
        dd_r = 0.5 * (r2 / dd_r + dd_r);
        dd_r = 0.5 * (r2 / dd_r + dd_r);
        double r = to_double(dd_r);

        r = ldexp_(r, exp);
        return r;
    }


    inline constexpr double sinh_(double x)
    {
        double a = expm1_(x);
        double b = expm1_(-x);
        return 0.5 * (a - b);
    }
    inline constexpr double cosh_(double x)
    {
        double a = expm1_(x);
        double b = expm1_(-x);
        return (1. + 0.5 * (a + b));
    }

    inline constexpr double atan_Taylor(double x)
    {
        // x - x^3/3 + x^5/5
        double t = x * x;
        double v_ = 0;

        constexpr int max_iters = 26;
        double a = max_iters * 2 + 3;
        double b = max_iters * 2 + 1;

        for (int i = 0; i < max_iters; ++i) {
            // v_ = v - 1
            v_ = (1. + v_) * a / b * t;
            a -= 2.;
            b -= 2.;
        }

        return x * (1 + v_);
    }

    inline constexpr double atan_(double x)
    {
        // tan(x - y) = (tanx - tany)/(1 + tanx tany)
        // x = tx + arctan((tx-ty£©/(1+txty))
        // arctan(x) = pi/2 - arctan(1/x)
        // 
        // arctan(0) = 0
        // arctan(inf) = pi/2

        bool sign = false;
        if (x < 0) {
            sign = true;
            x = -x;
        }

        constexpr double atan_talbe[] = {
            0,
            0.2449786631268641541720824812112758109141440983811840671273759146,
            0.4636476090008061162142562314612144020285370542861202638109330887,
            0.6435011087932843868028092287173226380415105911153123828656061187,
            0.7853981633974483096156608458198757210492923498437764552437361480,
        };

        constexpr double tan_table[] = {
            0, 0.25, 0.5, 0.75, 1.0,
        };

        bool inverse = false;
        if (x > 1.) {
            x = 1 / x;
            inverse = true;
        }

        double v = 0;
        double low = 1;
        size_t idx = 0;
        for (size_t i = 0; i < std::size(phi_table); ++i) {
            if (fb::abs_(phi_table[i] - x) < low) {
                low = fb::abs_(phi_table[i] - x);
                idx = i;
            }
        }
        double t = tan_table[idx];
        double ratio = (x - t) / (1 + t * x);

        double arctan_ratio = atan_Taylor(ratio);

        double r = atan_talbe[idx] + arctan_ratio;


        if (sign) return -r;
        else return r;
    }

    inline constexpr double atan(double x)
    {
        if (std::is_constant_evaluated()) {
            return atan_(x);
        } else {
            return std::atan(x);
        }
    }

}
#endif
#include <vector>
#include <random>

dd_real ddrand(void);
constexpr dd_real sqrt(const dd_real& a);

constexpr dd_real polyeval(const dd_real* c, int n, const dd_real& x);
constexpr dd_real polyroot(const dd_real* c, int n,
    const dd_real& x0, int max_iter = 32, double thresh = 0.0);

constexpr bool isnan(const dd_real& a) { return a.isnan(); }
constexpr bool isfinite(const dd_real& a) { return a.isfinite(); }
constexpr bool isinf(const dd_real& a) { return a.isinf(); }

QD_API constexpr dd_real rem(const dd_real& a, const dd_real& b);
QD_API constexpr dd_real drem(const dd_real& a, const dd_real& b);
QD_API constexpr dd_real divrem(const dd_real& a, const dd_real& b, dd_real& r);

QD_API constexpr dd_real pow(const dd_real& a, int n);
QD_API constexpr dd_real pow(const dd_real& a, const dd_real& b);

QD_API constexpr dd_real sqrt(const dd_real& a);
QD_API constexpr dd_real nroot(const dd_real& a, int n);

constexpr dd_real nint(const dd_real& a);
constexpr dd_real round(const dd_real& a);
constexpr dd_real floor(const dd_real& a);
constexpr dd_real ceil(const dd_real& a);
constexpr dd_real aint(const dd_real& a);

dd_real ddrand(void);

constexpr dd_real exp(const dd_real& a);
constexpr dd_real ldexp(const dd_real& a, int exp);
constexpr dd_real expm1(const dd_real& a);
constexpr dd_real log(const dd_real& a);
constexpr dd_real log10(const dd_real& a);

constexpr dd_real sin(const dd_real& a);
constexpr dd_real cos(const dd_real& a);
constexpr dd_real tan(const dd_real& a);
constexpr void sincos(const dd_real& a, dd_real& sin_a, dd_real& cos_a);

constexpr dd_real asin(const dd_real& a);
constexpr dd_real acos(const dd_real& a);
constexpr dd_real atan(const dd_real& a);
constexpr dd_real atan2(const dd_real& y, const dd_real& x);

constexpr dd_real sinh(const dd_real& a);
constexpr dd_real cosh(const dd_real& a);
constexpr dd_real tanh(const dd_real& a);
constexpr void sincosh(const dd_real& a,
    dd_real& sinh_a, dd_real& cosh_a);

constexpr dd_real asinh(const dd_real& a);
constexpr dd_real acosh(const dd_real& a);
constexpr dd_real atanh(const dd_real& a);

constexpr dd_real fabs(const dd_real& a);
constexpr dd_real abs(const dd_real& a);   /* same as fabs */

constexpr dd_real fmod(const dd_real& a, const dd_real& b);




/* Computes the square root of the double-double number dd.
   NOTE: dd must be a non-negative number.                   */
inline QD_CONSTEXPR dd_real sqrt(const dd_real &a) {
  /* Strategy:  Use Karp's trick:  if x is an approximation
     to sqrt(a), then

        sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)

     The approximation is accurate to twice the accuracy of x.
     Also, the multiplication (a*x) and [-]*x can be done with
     only half the precision.
  */

  if (a.is_zero())
    return a;

  if (a.is_negative()) {
    //dd_real::error("(dd_real::sqrt): Negative argument.");
    return dd_real::_nan;
  }

  double x = 1.0 / fb::sqrt(a.x[0]);
  double ax = a.x[0] * x;
  return dd_real::add(ax, (a - dd_real::sqr(ax)).x[0] * (x * 0.5));
}

/* Computes the square root of a double in double-double precision. 
   NOTE: d must not be negative.                                   */
inline QD_CONSTEXPR dd_real dd_real::sqrt(double d) {
  return ::sqrt(dd_real(d));
}

/* Computes the n-th root of the double-double number a.
   NOTE: n must be a positive integer.  
   NOTE: If n is even, then a must not be negative.       */
inline QD_CONSTEXPR dd_real nroot(const dd_real &a, int n) {
  /* Strategy:  Use Newton iteration for the function

          f(x) = x^(-n) - a

     to find its root a^{-1/n}.  The iteration is thus

          x' = x + x * (1 - a * x^n) / n

     which converges quadratically.  We can then find 
    a^{1/n} by taking the reciprocal.
  */

  if (n <= 0) {
    //dd_real::error("(dd_real::nroot): N must be positive.");
    return dd_real::_nan;
  }

  if (n%2 == 0 && a.is_negative()) {
    //dd_real::error("(dd_real::nroot): Negative argument.");
    return dd_real::_nan;
  }

  if (n == 1) {
    return a;
  } 
  if (n == 2) {
    return sqrt(a);
  }

  if (a.is_zero())
    return 0.0;

  /* Note  a^{-1/n} = exp(-log(a)/n) */
  dd_real r = abs(a);
  dd_real x = std::exp(-std::log(r.x[0]) / n);

  /* Perform Newton's iteration. */
  x += x * (1.0 - r * npwr(x, n)) / static_cast<double>(n);
  if (a.x[0] < 0.0)
    x = -x;
  return 1.0/x;
}

/********** Remainder **********/
inline QD_CONSTEXPR dd_real drem(const dd_real& a, const dd_real& b)
{
    dd_real n = nint(a / b);
    return (a - n * b);
}

inline QD_CONSTEXPR dd_real divrem(const dd_real& a, const dd_real& b, dd_real& r)
{
    dd_real n = nint(a / b);
    r = a - n * b;
    return n;
}


inline QD_CONSTEXPR dd_real pow(const dd_real &a, int n) {
  return npwr(a, n);
}

inline QD_CONSTEXPR dd_real pow(const dd_real &a, const dd_real &b) {
  return exp(b * log(a));
}

inline QD_CONSTEXPR int n_inv_fact_dd = 15;
inline QD_CONSTEXPR double inv_fact_dd[n_inv_fact_dd][2] = {
  { 1.66666666666666657e-01,  9.25185853854297066e-18},
  { 4.16666666666666644e-02,  2.31296463463574266e-18},
  { 8.33333333333333322e-03,  1.15648231731787138e-19},
  { 1.38888888888888894e-03, -5.30054395437357706e-20},
  { 1.98412698412698413e-04,  1.72095582934207053e-22},
  { 2.48015873015873016e-05,  2.15119478667758816e-23},
  { 2.75573192239858925e-06, -1.85839327404647208e-22},
  { 2.75573192239858883e-07,  2.37677146222502973e-23},
  { 2.50521083854417202e-08, -1.44881407093591197e-24},
  { 2.08767569878681002e-09, -1.20734505911325997e-25},
  { 1.60590438368216133e-10,  1.25852945887520981e-26},
  { 1.14707455977297245e-11,  2.06555127528307454e-28},
  { 7.64716373181981641e-13,  7.03872877733453001e-30},
  { 4.77947733238738525e-14,  4.39920548583408126e-31},
  { 2.81145725434552060e-15,  1.65088427308614326e-31}
};

/* Exponential.  Computes exp(x) in double-double precision. */
inline QD_CONSTEXPR dd_real exp(const dd_real &a) {
  /* Strategy:  We first reduce the size of x by noting that
     
          exp(kr + m * log(2)) = 2^m * exp(r)^k

     where m and k are integers.  By choosing m appropriately
     we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is 
     evaluated using the familiar Taylor series.  Reducing the 
     argument substantially speeds up the convergence.       */  

  const double k = 512.0;
  const double inv_k = 1.0 / k;

  if (a.x[0] <= -709.0)
    return 0.0;

  if (a.x[0] >=  709.0)
    return dd_real::_inf;

  if (a.is_zero())
    return 1.0;

  if (a.is_one())
    return dd_real::_e;

  double m = fb::round(a.x[0] / dd_real::_log2.x[0]);
  dd_real r = mul_pwr2(a - dd_real::_log2 * m, inv_k);
  dd_real s, t, p;

  p = sqr(r);
  s = r + mul_pwr2(p, 0.5);
  p *= r;
  t = p * dd_real(inv_fact_dd[0][0], inv_fact_dd[0][1]);
  int i = 0;
  do {
    s += t;
    p *= r;
    ++i;
    t = p * dd_real(inv_fact_dd[i][0], inv_fact_dd[i][1]);
  } while (fb::abs(to_double(t)) > inv_k * dd_real::_eps && i < 5);

  s += t;

  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s); 
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s += 1.0;

  return ldexp(s, static_cast<int>(m));
}

inline QD_CONSTEXPR dd_real ldexp(const dd_real& a, int exp)
{
    return dd_real(fb::ldexp(a.x[0], exp), fb::ldexp(a.x[1], exp));
}

inline QD_CONSTEXPR dd_real expm1(const dd_real& a)
{
    /* Strategy:  We first reduce the size of x by noting that

            exp(kr + m * log(2)) = 2^m * exp(r)^k

       where m and k are integers.  By choosing m appropriately
       we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
       evaluated using the familiar Taylor series.  Reducing the
       argument substantially speeds up the convergence.       */

    const double k = 512.0;
    const double inv_k = 1.0 / k;

    if (a.x[0] <= -709.0)
        return 0.0;

    if (a.x[0] >= 709.0)
        return dd_real::_inf;

    if (a.is_zero())
        return 1.0;

    if (a.is_one())
        return dd_real::_e;

    double m = fb::round(a.x[0] / dd_real::_log2.x[0]);
    dd_real r = mul_pwr2(a - dd_real::_log2 * m, inv_k);
    dd_real s, t, p;

    p = sqr(r);
    s = r + mul_pwr2(p, 0.5);
    p *= r;
    t = p * dd_real(inv_fact_dd[0][0], inv_fact_dd[0][1]);
    int i = 0;
    do {
        s += t;
        p *= r;
        ++i;
        t = p * dd_real(inv_fact_dd[i][0], inv_fact_dd[i][1]);
    } while (fb::abs(to_double(t)) > inv_k * dd_real::_eps && i < 5);

    s += t;

    // (1+s)^2
    // (1+s)^2 = s^2 + 2s + 1 =  1 + s'
    // (1+s)^2 = 1 + s'
    // s' = s^2 + 2s
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);
    s = mul_pwr2(s, 2.0) + sqr(s);

    // (1+s)2^m - 1
    if (m == 0) {
        return s;
    } else if (m == 1) {
        return 1. + mul_pwr2(s, 2.);
    } else {
        return ldexp(1. + s, static_cast<int>(m)) - 1;
    }
}

inline QD_CONSTEXPR dd_real log1p_smallfrac(const dd_real& frac)
{
    dd_real x = fb::log1p(frac.x[0]);   /* Initial approximation */
    dd_real expm1_ = expm1(x);
    dd_real a_1 = frac;
    x = x + (a_1 - expm1_) / (1. + expm1_);
    return x;
}

/* Logarithm.  Computes log(x) in double-double precision.
   This is a natural logarithm (i.e., base e).            */
inline QD_CONSTEXPR dd_real log(const dd_real &a) {
  /* Strategy.  The Taylor series for log converges much more
     slowly than that of exp, due to the lack of the factorial
     term in the denominator.  Hence this routine instead tries
     to determine the root of the function

         f(x) = exp(x) - a

     using Newton iteration.  The iteration is given by

         x' = x - f(x)/f'(x) 
            = x - (1 - a * exp(-x))
            = x + a * exp(-x) - 1.
           
     Only one iteration is needed, since Newton's iteration
     approximately doubles the number of digits per iteration. */

  if (a.is_one()) {
    return 0.0;
  }

  if (a.x[0] <= 0.0) {
    return dd_real::_nan;
  }


  if (fb::abs(a.x[0] - 1.) < 0.5) {
      return log1p_smallfrac(a - 1.);
  } else {
      dd_real x = fb::log(a.x[0]);   /* Initial approximation */
      x = x + a * exp(-x) - 1.0;
      return x;
  }

}

inline QD_CONSTEXPR dd_real log1p(const dd_real& a)
{
    dd_real aabs = fabs(a);
    if (aabs < 0.5) {
        return log1p_smallfrac(a);
    } else {
        return log(1. + a);
    }
}

inline QD_CONSTEXPR dd_real log10(const dd_real &a) {
  return log(a) / dd_real::_log10;
}

inline QD_CONSTEXPR const dd_real _pi16 = dd_real(1.963495408493620697e-01,
                                     7.654042494670957545e-18);

/* Table of sin(k * pi/16) and cos(k * pi/16). */
inline QD_CONSTEXPR const double sin_table_dd [4][2] = {
  {1.950903220161282758e-01, -7.991079068461731263e-18},
  {3.826834323650897818e-01, -1.005077269646158761e-17},
  {5.555702330196021776e-01,  4.709410940561676821e-17},
  {7.071067811865475727e-01, -4.833646656726456726e-17}
};

inline QD_CONSTEXPR const double cos_table_dd [4][2] = {
  {9.807852804032304306e-01, 1.854693999782500573e-17},
  {9.238795325112867385e-01, 1.764504708433667706e-17},
  {8.314696123025452357e-01, 1.407385698472802389e-18},
  {7.071067811865475727e-01, -4.833646656726456726e-17}
};

/* Computes sin(a) using Taylor series.
   Assumes |a| <= pi/32.                           */
inline QD_CONSTEXPR dd_real sin_taylor(const dd_real &a) {
  const double thresh = 0.5 * fb::abs(to_double(a)) * dd_real::_eps;
  dd_real r, s, t, x;

  if (a.is_zero()) {
    return 0.0;
  }

  int i = 0;
  x = -sqr(a);
  s = a;
  r = a;
  do {
    r *= x;
    t = r * dd_real(inv_fact_dd[i][0], inv_fact_dd[i][1]);
    s += t;
    i += 2;
  } while (i < n_inv_fact_dd && fb::abs(to_double(t)) > thresh);

  return s;
}

inline QD_CONSTEXPR dd_real cos_taylor(const dd_real &a) {
  const double thresh = 0.5 * dd_real::_eps;
  dd_real r, s, t, x;

  if (a.is_zero()) {
    return 1.0;
  }

  x = -sqr(a);
  r = x;
  s = 1.0 + mul_pwr2(r, 0.5);
  int i = 1;
  do {
    r *= x;
    t = r * dd_real(inv_fact_dd[i][0], inv_fact_dd[i][1]);
    s += t;
    i += 2;
  } while (i < n_inv_fact_dd && fb::abs(to_double(t)) > thresh);

  return s;
}

inline QD_CONSTEXPR void sincos_taylor(const dd_real &a,
                          dd_real &sin_a, dd_real &cos_a) {
  if (a.is_zero()) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  sin_a = sin_taylor(a);
  cos_a = sqrt(1.0 - sqr(sin_a));
}


inline QD_CONSTEXPR dd_real sin(const dd_real &a) {

  /* Strategy.  To compute sin(x), we choose integers a, b so that

       x = s + a * (pi/2) + b * (pi/16)

     and |s| <= pi/32.  Using the fact that 

       sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))

     we can compute sin(x) from sin(s), cos(s).  This greatly 
     increases the convergence of the sine Taylor series. */

  if (a.is_zero()) {
    return 0.0;
  }

  // approximately reduce modulo 2*pi
  dd_real z = nint(a / dd_real::_2pi);
  dd_real r = a - dd_real::_2pi * z;

  // approximately reduce modulo pi/2 and then modulo pi/16.
  dd_real t;
  double q = fb::round(r.x[0] / dd_real::_pi2.x[0]);
  t = r - dd_real::_pi2 * q;
  int j = static_cast<int>(q);
  q = fb::round(t.x[0] / _pi16.x[0]);
  t -= _pi16 * q;
  int k = static_cast<int>(q);
  int abs_k = fb::abs(k);

  if (j < -2 || j > 2) {
    //dd_real::error("(dd_real::sin): Cannot reduce modulo pi/2.");
    return dd_real::_nan;
  }

  if (abs_k > 4) {
    //dd_real::error("(dd_real::sin): Cannot reduce modulo pi/16.");
    return dd_real::_nan;
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return sin_taylor(t);
      case 1:
        return cos_taylor(t);
      case -1:
        return -cos_taylor(t);
      default:
        return -sin_taylor(t);
    }
  }

  dd_real u(cos_table_dd[abs_k-1][0], cos_table_dd[abs_k-1][1]);
  dd_real v(sin_table_dd[abs_k-1][0], sin_table_dd[abs_k-1][1]);
  dd_real sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);
  if (j == 0) {
    if (k > 0) {
      r = u * sin_t + v * cos_t;
    } else {
      r = u * sin_t - v * cos_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      r = u * cos_t - v * sin_t;
    } else {
      r = u * cos_t + v * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      r = v * sin_t - u * cos_t;
    } else if (k < 0) {
      r = -u * cos_t - v * sin_t;
    }
  } else {
    if (k > 0) {
      r = -u * sin_t - v * cos_t;
    } else {
      r = v * cos_t - u * sin_t;
    }
  }

  return r;
}

inline QD_CONSTEXPR dd_real cos(const dd_real &a) {

  if (a.is_zero()) {
    return 1.0;
  }

  // approximately reduce modulo 2*pi
  dd_real z = nint(a / dd_real::_2pi);
  dd_real r = a - z * dd_real::_2pi;

  // approximately reduce modulo pi/2 and then modulo pi/16
  dd_real t;
  double q = fb::round(r.x[0] / dd_real::_pi2.x[0]);
  t = r - dd_real::_pi2 * q;
  int j = static_cast<int>(q);
  q = fb::round(t.x[0] / _pi16.x[0]);
  t -= _pi16 * q;
  int k = static_cast<int>(q);
  int abs_k = fb::abs(k);

  if (j < -2 || j > 2) {
    //dd_real::error("(dd_real::cos): Cannot reduce modulo pi/2.");
    return dd_real::_nan;
  }

  if (abs_k > 4) {
    //dd_real::error("(dd_real::cos): Cannot reduce modulo pi/16.");
    return dd_real::_nan;
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return cos_taylor(t);
      case 1:
        return -sin_taylor(t);
      case -1:
        return sin_taylor(t);
      default:
        return -cos_taylor(t);
    }
  }

  dd_real sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);
  dd_real u(cos_table_dd[abs_k-1][0], cos_table_dd[abs_k-1][1]);
  dd_real v(sin_table_dd[abs_k-1][0], sin_table_dd[abs_k-1][1]);

  if (j == 0) {
    if (k > 0) {
      r = u * cos_t - v * sin_t;
    } else {
      r = u * cos_t + v * sin_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      r = - u * sin_t - v * cos_t;
    } else {
      r = v * cos_t - u * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      r = u * sin_t + v * cos_t;
    } else {
      r = u * sin_t - v * cos_t;
    }
  } else {
    if (k > 0) {
      r = v * sin_t - u * cos_t;
    } else {
      r = - u * cos_t - v * sin_t;
    }
  }

  return r;
}

inline QD_CONSTEXPR void sincos(const dd_real &a, dd_real &sin_a, dd_real &cos_a) {

  if (a.is_zero()) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  // approximately reduce modulo 2*pi
  dd_real z = nint(a / dd_real::_2pi);
  dd_real r = a - dd_real::_2pi * z;

  // approximately reduce module pi/2 and pi/16
  dd_real t;
  double q = fb::floor(r.x[0] / dd_real::_pi2.x[0] + 0.5);
  t = r - dd_real::_pi2 * q;
  int j = static_cast<int>(q);
  int abs_j = fb::abs(j);
  q = fb::floor(t.x[0] / _pi16.x[0] + 0.5);
  t -= _pi16 * q;
  int k = static_cast<int>(q);
  int abs_k = fb::abs(k);

  if (abs_j > 2) {
    //dd_real::error("(dd_real::sincos): Cannot reduce modulo pi/2.");
    cos_a = sin_a = dd_real::_nan;
    return;
  }

  if (abs_k > 4) {
    //dd_real::error("(dd_real::sincos): Cannot reduce modulo pi/16.");
    cos_a = sin_a = dd_real::_nan;
    return;
  }

  dd_real sin_t, cos_t;
  dd_real s, c;

  sincos_taylor(t, sin_t, cos_t);

  if (abs_k == 0) {
    s = sin_t;
    c = cos_t;
  } else {
    dd_real u(cos_table_dd[abs_k-1][0], cos_table_dd[abs_k-1][1]);
    dd_real v(sin_table_dd[abs_k-1][0], sin_table_dd[abs_k-1][1]);

    if (k > 0) {
      s = u * sin_t + v * cos_t;
      c = u * cos_t - v * sin_t;
    } else {
      s = u * sin_t - v * cos_t;
      c = u * cos_t + v * sin_t;
    }
  }

  if (abs_j == 0) {
    sin_a = s;
    cos_a = c;
  } else if (j == 1) {
    sin_a = c;
    cos_a = -s;
  } else if (j == -1) {
    sin_a = -c;
    cos_a = s;
  } else {
    sin_a = -s;
    cos_a = -c;
  }
  
}

inline QD_CONSTEXPR dd_real atan(const dd_real &a) {
  return atan2(a, dd_real(1.0));
}

inline dd_real QD_CONSTEXPR atan2(const dd_real &y, const dd_real &x) {
  /* Strategy: Instead of using Taylor series to compute 
     arctan, we instead use Newton's iteration to solve
     the equation

        sin(z) = y/r    or    cos(z) = x/r

     where r = sqrt(x^2 + y^2).
     The iteration is given by

        z' = z + (y - sin(z)) / cos(z)          (for equation 1)
        z' = z - (x - cos(z)) / sin(z)          (for equation 2)

     Here, x and y are normalized so that x^2 + y^2 = 1.
     If |x| > |y|, then first iteration is used since the 
     denominator is larger.  Otherwise, the second is used.
  */

  if (x.is_zero()) {
    
    if (y.is_zero()) {
      /* Both x and y is zero. */
      //dd_real::error("(dd_real::atan2): Both arguments zero.");
      return dd_real::_nan;
    }

    return (y.is_positive()) ? dd_real::_pi2 : -dd_real::_pi2;
  } else if (y.is_zero()) {
    return (x.is_positive()) ? dd_real(0.0) : dd_real::_pi;
  }

  if (x == y) {
    return (y.is_positive()) ? dd_real::_pi4 : -dd_real::_3pi4;
  }

  if (x == -y) {
    return (y.is_positive()) ? dd_real::_3pi4 : -dd_real::_pi4;
  }

  dd_real r = sqrt(sqr(x) + sqr(y));
  dd_real xx = x / r;
  dd_real yy = y / r;

  /* Compute double precision approximation to atan. */
  dd_real z = std::atan2(to_double(y), to_double(x));
  dd_real sin_z, cos_z;

  if (fb::abs(xx.x[0]) > fb::abs(yy.x[0])) {
    /* Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z)  */
    sincos(z, sin_z, cos_z);
    z += (yy - sin_z) / cos_z;
  } else {
    /* Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z)  */
    sincos(z, sin_z, cos_z);
    z -= (xx - cos_z) / sin_z;
  }

  return z;
}

inline QD_CONSTEXPR dd_real tan(const dd_real &a) {
  dd_real s, c;
  sincos(a, s, c);
  return s/c;
}

inline QD_CONSTEXPR dd_real asin(const dd_real &a) {
  dd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    return dd_real::_nan;
  }

  if (abs_a.is_one()) {
    return (a.is_positive()) ? dd_real::_pi2 : -dd_real::_pi2;
  }

  return atan2(a, sqrt(1.0 - sqr(a)));
}

inline QD_CONSTEXPR dd_real acos(const dd_real &a) {
  dd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    //dd_real::error("(dd_real::acos): Argument out of domain.");
    return dd_real::_nan;
  }

  if (abs_a.is_one()) {
    return (a.is_positive()) ? dd_real(0.0) : dd_real::_pi;
  }

  return atan2(sqrt(1.0 - sqr(a)), a);
}
 
inline QD_CONSTEXPR dd_real sinh(const dd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  if (abs(a) > 0.05) {
    dd_real ea = exp(a);
    return mul_pwr2(ea - inv(ea), 0.5);
  }

  /* since a is small, using the above formula gives
     a lot of cancellation.  So use Taylor series.   */
  dd_real s = a;
  dd_real t = a;
  dd_real r = sqr(t);
  double m = 1.0;
  double thresh = fb::abs((to_double(a)) * dd_real::_eps);

  do {
    m += 2.0;
    t *= r;
    t /= (m-1) * m;

    s += t;
  } while (abs(t) > thresh);

  return s;

}

inline dd_real QD_CONSTEXPR cosh(const dd_real &a) {
  if (a.is_zero()) {
    return 1.0;
  }

  dd_real ea = exp(a);
  return mul_pwr2(ea + inv(ea), 0.5);
}

dd_real QD_CONSTEXPR tanh(const dd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  dd_real ea = exp(a);
  dd_real inv_ea = inv(ea);
  return (ea - inv_ea) / (ea + inv_ea);
}

inline void QD_CONSTEXPR sincosh(const dd_real &a, dd_real &sinh_a, dd_real &cosh_a) {
  sinh_a = sinh(a);
  cosh_a = cosh(a);
}

inline dd_real QD_CONSTEXPR asinh(const dd_real &a) {
    dd_real a2 = sqr(a);
    dd_real r = a2 + 1.0;
    if (a2 <= 0.25) {
        if (a >= 0.)
            return log1p(a + a2 / (1. + sqrt(r)));
        else
            return -log1p(-a + a2 / (1. + sqrt(r)));
    } else {
        if (a >= 0.)
            return log(a + sqrt(r));
        else
            return -log(-a + sqrt(r));
    }
}

inline dd_real QD_CONSTEXPR acosh(const dd_real &a) {
    if (a < 1.0) {
        return dd_real::_nan;
    }
    return log(a + sqrt(sqr(a) - 1.0));
}

inline QD_CONSTEXPR dd_real atanh(const dd_real &a) {
    if (abs(a) >= 1.0) {
        return dd_real::_nan;
    }
    return mul_pwr2(log((1.0 + a) / (1.0 - a)), 0.5);
}

inline QD_CONSTEXPR dd_real fmod(const dd_real &a, const dd_real &b) {
  dd_real n = aint(a / b);
  return (a - b * n);
}


/* Round to Nearest integer */
inline QD_CONSTEXPR dd_real nint(const dd_real& a)
{
    double hi = fb::round(a.x[0]);
    double lo;

    if (hi == a.x[0]) {
        /* High word is an integer already.  Round the low word.*/
        lo = fb::round(a.x[1]);

        /* Renormalize. This is needed if x[0] = some integer, x[1] = 1/2.*/
        hi = qd::quick_two_sum(hi, lo, lo);
    } else {
        /* High word is not an integer. */
        lo = 0.0;
        if (fb::abs(hi - a.x[0]) == 0.5 && a.x[1] < 0.0) {
            /* There is a tie in the high word, consult the low word
               to break the tie. */
            hi -= 1.0;      /* NOTE: This does not cause INEXACT. */
        }
    }

    return dd_real(hi, lo);
}

inline QD_CONSTEXPR dd_real round(const dd_real& a)
{
    return nint(a);
}

inline constexpr dd_real floor(const dd_real& a)
{
    double hi = fb::floor(a.x[0]);
    double lo = 0.0;

    if (hi == a.x[0]) {
        /* High word is integer already.  Round the low word. */
        lo = fb::floor(a.x[1]);
        hi = qd::quick_two_sum(hi, lo, lo);
    }

    return dd_real(hi, lo);
}

inline constexpr dd_real ceil(const dd_real& a)
{
    double hi = fb::ceil(a.x[0]);
    double lo = 0.0;

    if (hi == a.x[0]) {
        /* High word is integer already.  Round the low word. */
        lo = std::ceil(a.x[1]);
        hi = qd::quick_two_sum(hi, lo, lo);
    }

    return dd_real(hi, lo);
}

inline constexpr dd_real aint(const dd_real& a)
{
    return (a.x[0] >= 0.0) ? floor(a) : ceil(a);
}

/* Random number generator */
inline dd_real dd_real::rand()
{
    return ddrand();
}


inline std::mt19937 mt19937;
inline dd_real ddrand() {
  static const double m_const = 4.6566128730773926e-10;  /* = 2^{-31} */
  double m = m_const;
  dd_real r = 0.0;
  double d;

  /* Strategy:  Generate 31 bits at a time, using lrand48 
     random number generator.  Shift the bits, and reapeat
     4 times. */
  // lrand48()   not partable
  // std::rand() too bad quality

  for (int i = 0; i < 4; i++, m *= m_const) {
    d = mt19937() * m;
    r += d;
  }

  return r;
}

/* polyeval(c, n, x)
   Evaluates the given n-th degree polynomial at x.
   The polynomial is given by the array of (n+1) coefficients. */
inline QD_CONSTEXPR dd_real polyeval(const dd_real *c, int n, const dd_real &x) {
  /* Just use Horner's method of polynomial evaluation. */
  dd_real r = c[n];
  
  for (int i = n-1; i >= 0; i--) {
    r *= x;
    r += c[i];
  }

  return r;
}

/* polyroot(c, n, x0)
   Given an n-th degree polynomial, finds a root close to 
   the given guess x0.  Note that this uses simple Newton
   iteration scheme, and does not work for multiple roots.  */
inline QD_CONSTEXPR dd_real polyroot(const dd_real *c, int n,
    const dd_real &x0, int max_iter, double thresh) {
  dd_real x = x0;
  dd_real f;
  dd_real *d = new dd_real[n]; //HL: safe new
  bool conv = false;
  int i;
  double max_c = fb::abs(to_double(c[0]));
  double v;

  if (thresh == 0.0) thresh = dd_real::_eps;

  /* Compute the coefficients of the derivatives. */
  for (i = 1; i <= n; i++) {
    v = std::abs(to_double(c[i]));
    if (v > max_c) max_c = v;
    d[i-1] = c[i] * static_cast<double>(i);
  }
  thresh *= max_c;

  /* Newton iteration. */
  for (i = 0; i < max_iter; i++) {
    f = polyeval(c, n, x);

    if (abs(f) < thresh) {
      conv = true;
      break;
    }
    x -= (f / polyeval(d, n-1, x));
  }
  delete [] d;

  if (!conv) {
    //dd_real::error("(dd_real::polyroot): Failed to converge.");
    return dd_real::_nan;
  }

  return x;
}

#if 0
/* Constructor.  Reads a double-double number from the string s
   and constructs a double-double number.                         */
dd_real::dd_real(const char *s) {
  if (dd_real::read(s, *this)) {
    dd_real::error("(dd_real::dd_real): INPUT ERROR.");
    *this = dd_real::_nan;
  }
}

dd_real &dd_real::operator=(const char *s) {
  if (dd_real::read(s, *this)) {
    dd_real::error("(dd_real::operator=): INPUT ERROR.");
    *this = dd_real::_nan;
  }
  return *this;
}



inline dd_real dd_real::debug_rand()
{

    if (std::rand() % 2 == 0)
        return ddrand();

    int expn = 0;
    dd_real a = 0.0;
    double d;
    for (int i = 0; i < 2; i++) {
        d = std::ldexp(static_cast<double>(std::rand()) / RAND_MAX, -expn);
        a += d;
        expn = expn + 54 + std::rand() % 200;
    }
    return a;
}

#endif

#endif

/*  dd_io.h  */
#ifndef _QD_DD_IO_H
#define _QD_DD_IO_H


std::ostream& operator<<(std::ostream& s, const dd_real& a);
std::istream& operator>>(std::istream& s, dd_real& a);

inline dd_real ldexp_(const dd_real& a, int exp)
{
    return dd_real(std::ldexp(a.x[0], exp), std::ldexp(a.x[1], exp));
}

/* Outputs the double-double number dd. */
inline std::ostream& operator<<(std::ostream& os, const dd_real& dd)
{
    bool showpos = (os.flags() & std::ios_base::showpos) != 0;
    bool uppercase = (os.flags() & std::ios_base::uppercase) != 0;
    return os << dd.to_string((int)/*HL: */os.precision(), (int)/*HL: */os.width(), os.flags(),
        showpos, uppercase, os.fill());
}

/* Reads in the double-double number a. */
inline std::istream& operator>>(std::istream& s, dd_real& a)
{
    char str[255];
    s >> str;
    dd_real::read(str, a);
    return s;
}

inline void dd_real::to_digits(char* s, int& expn, int precision) const
{
    int D = precision + 1;  /* number of digits to compute */

    dd_real r = abs(*this);
    int e;  /* exponent */
    int i, d;

    if (x[0] == 0.0) {
        /* this == 0.0 */
        expn = 0;
        for (i = 0; i < precision; i++) s[i] = '0';
        return;
    }

    /* First determine the (approximate) exponent. */
    e = to_int(std::floor(std::log10(std::abs(x[0]))));

    if (e < -300) {
        r *= dd_real(10.0) ^ 300;
        r /= dd_real(10.0) ^ (e + 300);
    } else if (e > 300) {
        r = ldexp_(r, -53);
        r /= dd_real(10.0) ^ e;
        r = ldexp_(r, 53);
    } else {
        r /= dd_real(10.0) ^ e;
    }

    /* Fix exponent if we are off by one */
    if (r >= 10.0) {
        r /= 10.0;
        e++;
    } else if (r < 1.0) {
        r *= 10.0;
        e--;
    }

    if (r >= 10.0 || r < 1.0) {
        qd::error("(dd_real::to_digits): can't compute exponent.");
        return;
    }

    /* Extract the digits */
    for (i = 0; i < D; i++) {
        d = static_cast<int>(r.x[0]);
        r -= d;
        r *= 10.0;

        s[i] = static_cast<char>(d + '0');
    }

    /* Fix out of range digits. */
    for (i = D - 1; i > 0; i--) {
        if (s[i] < '0') {
            s[i - 1]--;
            s[i] += 10;
        } else if (s[i] > '9') {
            s[i - 1]++;
            s[i] -= 10;
        }
    }

    if (s[0] <= '0') {
        qd::error("(dd_real::to_digits): non-positive leading digit.");
        return;
    }

    /* Round, handle carry */
    if (s[D - 1] >= '5') {
        s[D - 2]++;

        i = D - 2;
        while (i > 0 && s[i] > '9') {
            s[i] -= 10;
            s[--i]++;
        }
    }

    /* If first digit is 10, shift everything. */
    if (s[0] > '9') {
        e++;
        for (i = precision; i >= 2; i--) s[i] = s[i - 1];
        s[0] = '1';
        s[1] = '0';
    }

    s[precision] = 0;
    expn = e;
}

/* Writes the double-double number into the character array s of length len.
   The integer d specifies how many significant digits to write.
   The string s must be able to hold at least (d+8) characters.
   showpos indicates whether to use the + sign, and uppercase indicates
   whether the E or e is to be used for the exponent. */
inline void dd_real::write(char* s, int len, int precision,
    bool showpos, bool uppercase) const
{
    std::string str = to_string(precision, 0, std::ios_base::scientific, showpos, uppercase);

    // MSVC don't like user to use strncpy
    size_t ln = (std::min)(str.length(), (size_t)len - 1);
    memcpy(s, str.data(), ln);
    s[len - 1] = '\0';
}


inline void round_string(char* s, int precision, int* offset)
{
    /*
     Input string must be all digits or errors will occur.
     */

    int i;
    int D = precision;

    /* Round, handle carry */
    if (s[D - 1] >= '5') {
        s[D - 2]++;

        i = D - 2;
        while (i > 0 && s[i] > '9') {
            s[i] -= 10;
            s[--i]++;
        }
    }

    /* If first digit is 10, shift everything. */
    if (s[0] > '9') {
        // e++; // don't modify exponent here
        for (i = precision; i >= 2; i--) s[i] = s[i - 1];
        s[0] = '1';
        s[1] = '0';

        (*offset)++; // now offset needs to be increased by one
        precision++;
    }

    s[precision] = 0; // add terminator for array
}

inline std::string dd_real::to_string(int precision, int width, std::ios_base::fmtflags fmt,
    bool showpos, bool uppercase, char fill) const
{
    std::string s;
    bool fixed = (fmt & std::ios_base::fixed) != 0;
    bool sgn = true;
    int i, e = 0;

    if (isnan()) {
        s = uppercase ? "NAN" : "nan";
        sgn = false;
    } else {
        if (*this < 0.0)
            s += '-';
        else if (showpos)
            s += '+';
        else
            sgn = false;

        if (isinf()) {
            s += uppercase ? "INF" : "inf";
        } else if (*this == 0.0) {
            /* Zero case */
            s += '0';
            if (precision > 0) {
                s += '.';
                s.append(precision, '0');
            }
        } else {
            /* Non-zero case */
            int off = (fixed ? (1 + to_int(floor(log10(abs(*this))))) : 1);
            int d = precision + off;

            int d_with_extra = d;
            if (fixed)
                d_with_extra = (std::max)(60, d); // longer than the max accuracy for DD

            // highly special case - fixed mode, precision is zero, abs(*this) < 1.0
            // without this trap a number like 0.9 printed fixed with 0 precision prints as 0
            // should be rounded to 1.
            if (fixed && (precision == 0) && (abs(*this) < 1.0)) {
                if (abs(*this) >= 0.5)
                    s += '1';
                else
                    s += '0';

                return s;
            }

            // handle near zero to working precision (but not exactly zero)
            if (fixed && d <= 0) {
                s += '0';
                if (precision > 0) {
                    s += '.';
                    s.append(precision, '0');
                }
            } else { // default

                char* t; //  = new char[d+1];
                int j;

                if (fixed) {
                    t = new char[(size_t)d_with_extra + 1]; //HL: safe new
                    to_digits(t, e, d_with_extra);
                } else {
                    t = new char[(size_t)d + 1];   //HL: safe new
                    to_digits(t, e, d);
                }

                if (fixed) {
                    // fix the string if it's been computed incorrectly
                    // round here in the decimal string if required
                    round_string(t, d + 1, &off);

                    if (off > 0) {
                        for (i = 0; i < off; i++) s += t[i];
                        if (precision > 0) {
                            s += '.';
                            for (j = 0; j < precision; j++, i++) s += t[i];
                        }
                    } else {
                        s += "0.";
                        if (off < 0) s.append(-off, '0');
                        for (i = 0; i < d; i++) s += t[i];
                    }
                } else {
                    s += t[0];
                    if (precision > 0) s += '.';

                    for (i = 1; i <= precision; i++)
                        s += t[i];

                }
                delete[] t;
            }
        }

        // trap for improper offset with large values
        // without this trap, output of values of the for 10^j - 1 fail for j > 28
        // and are output with the point in the wrong place, leading to a dramatically off value
        if (fixed && (precision > 0)) {
            // make sure that the value isn't dramatically larger
            double from_string = atof(s.c_str());

            // if this ratio is large, then we've got problems
            if (fabs(from_string / this->x[0]) > 3.0) {

                //HL: unreferenced local variable   int point_position;
                //HL: unreferenced local variable   char temp;

                // loop on the string, find the point, move it up one
                // don't act on the first character
                for (i = 1; i < (int)/*HL: */s.length(); i++) {
                    if (s[i] == '.') {
                        s[i] = s[i - 1];
                        s[i - 1] = '.';
                        break;
                    }
                }

                from_string = atof(s.c_str());
                // if this ratio is large, then the string has not been fixed
                if (fabs(from_string / this->x[0]) > 3.0) {
                    qd::error("Re-rounding unsuccessful in large number fixed point trap.");
                }
            }
        }


        if (!fixed && !isinf()) {
            /* Fill in exponent part */
            s += uppercase ? 'E' : 'e';
            append_expn(s, e);
        }
    }

    /* Fill in the blanks */
    int len = (int)s.length(); //HL:suppress warnings
    if (len < width) {
        int delta = width - len;
        if (fmt & std::ios_base::internal) {
            if (sgn)
                s.insert(static_cast<std::string::size_type>(1), delta, fill);
            else
                s.insert(static_cast<std::string::size_type>(0), delta, fill);
        } else if (fmt & std::ios_base::left) {
            s.append(delta, fill);
        } else {
            s.insert(static_cast<std::string::size_type>(0), delta, fill);
        }
    }

    return s;
}


/* Debugging routines */
inline void dd_real::dump(const std::string& name, std::ostream& os) const
{
    std::ios_base::fmtflags old_flags = os.flags();
    std::streamsize old_prec = os.precision(19);
    os << std::scientific;

    if (name.length() > 0) os << name << " = ";
    os << "[ " << std::setw(27) << x[0] << ", " << std::setw(27) << x[1] << " ]" << std::endl;

    os.precision(old_prec);
    os.flags(old_flags);
}

inline void dd_real::dump_bits(const std::string& name, std::ostream& os) const
{
    std::string::size_type len = name.length();
    if (len > 0) {
        os << name << " = ";
        len += 3;
    }
    os << "[ ";
    len += 2;
    print_double_info(os, x[0]);
    os << std::endl;
    for (std::string::size_type i = 0; i < len; i++) os << ' ';
    print_double_info(os, x[1]);
    os << " ]" << std::endl;
}
#endif

#endif
/*  qd_real.h  */
/*
 * include/qd_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Quad-double precision (>= 212-bit significand) floating point arithmetic
 * package, written in ANSI C++, taking full advantage of operator overloading.
 * Uses similar techniques as that of David Bailey's double-double package 
 * and that of Jonathan Shewchuk's adaptive precision floating point 
 * arithmetic package.  See
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *   http://www.cs.cmu.edu/~quake/robust.html
 *
 * for more details.
 *
 * Yozo Hida
 */
#ifndef _QD_QD_REAL_H
#define _QD_QD_REAL_H

#include <iostream>
#include <string>
#include <limits>
/*  fpu.h  */
/*
 * include/fpu.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 * Contains functions to set and restore the round-to-double flag in the
 * control word of a x86 FPU.  The algorithms in the double-double and
 * quad-double package does not function with the extended mode found in
 * these FPU.
 */
#ifndef _QD_FPU_H
#define _QD_FPU_H

 /*
  * Set the round-to-double flag, and save the old control word in old_cw.
  * If old_cw is NULL, the old control word is not saved.
  */
void fpu_fix_start(unsigned int* old_cw);

/*
 * Restore the control word.
 */
void fpu_fix_end(unsigned int* old_cw);



#ifdef X86
#ifdef  _WIN32
#include <float.h>
#else

#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

#ifndef _FPU_GETCW
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

#ifndef _FPU_DOUBLE
#define _FPU_DOUBLE 0x0200
#endif

#endif
#endif /* X86 */


inline void fpu_fix_start(unsigned int* old_cw)
{
#ifdef X86
#ifdef _WIN32
#ifndef _M_X64 

    /* Win 32 MSVC */
    unsigned int cw = _control87(0, 0);
    _control87(0x00010000, 0x00030000);
    if (old_cw) {
        *old_cw = cw;
    }

#endif
#else
    /* Linux */
    volatile unsigned short cw, new_cw;
    _FPU_GETCW(cw);

    new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(new_cw);

    if (old_cw) {
        *old_cw = cw;
    }
#endif
#endif
}

inline void fpu_fix_end(unsigned int* old_cw)
{
#ifdef X86
#ifdef _WIN32
#ifndef _M_X64 

    /* Win 32 MSVC */
    if (old_cw) {
        _control87(*old_cw, 0xFFFFFFFF);
    }

#endif // #ifndef _M_X64 

#else
    /* Linux */
    if (old_cw) {
        int cw;
        cw = *old_cw;
        _FPU_SETCW(cw);
    }
#endif
#endif
}

#endif  /* _QD_FPU_H */

struct QD_API qd_real {
  double x[4];    /* The Components. */

  /* Eliminates any zeros in the middle component(s). */
  void zero_elim();
  void zero_elim(double &e);

  QD_CONSTEXPR void renorm();
  QD_CONSTEXPR void renorm(double &e);

  void quick_accum(double d, double &e);
  void quick_prod_accum(double a, double b, double &e);

  QD_CONSTEXPR qd_real(double x0, double x1, double x2, double x3);
  explicit qd_real(const double *xx);

  static const qd_real _2pi;
  static const qd_real _pi;
  static const qd_real _3pi4;
  static const qd_real _pi2;
  static const qd_real _pi4;
  static const qd_real _e;
  static const qd_real _log2;
  static const qd_real _log10;
  static const qd_real _nan;
  static const qd_real _inf;

  static const double _eps;
  static const double _min_normalized;
  static const qd_real _max;
  static const qd_real _safe_max;
  static const int _ndigits;

  QD_CONSTEXPR qd_real();
#if 0
  qd_real(const char *s);
#endif
  explicit QD_CONSTEXPR qd_real(const dd_real &dd);
  QD_CONSTEXPR qd_real(double d);
  QD_CONSTEXPR qd_real(int i);

  QD_CONSTEXPR double operator[](int i) const;
  QD_CONSTEXPR double &operator[](int i);

  static void error(const char *msg);

  QD_CONSTEXPR bool isnan() const;
  QD_CONSTEXPR bool isfinite() const { return QD_ISFINITE(x[0]); }
  QD_CONSTEXPR bool isinf() const { return QD_ISINF(x[0]); }

  QD_CONSTEXPR static qd_real ieee_add(const qd_real &a, const qd_real &b);
  QD_CONSTEXPR static qd_real sloppy_add(const qd_real &a, const qd_real &b);

  QD_CONSTEXPR qd_real &operator+=(double a);
  QD_CONSTEXPR qd_real &operator+=(const dd_real &a);
  QD_CONSTEXPR qd_real &operator+=(const qd_real &a);

  QD_CONSTEXPR qd_real &operator-=(double a);
  QD_CONSTEXPR qd_real &operator-=(const dd_real &a);
  QD_CONSTEXPR qd_real &operator-=(const qd_real &a);

  QD_CONSTEXPR static qd_real sloppy_mul(const qd_real &a, const qd_real &b);
  QD_CONSTEXPR static qd_real accurate_mul(const qd_real &a, const qd_real &b);

  QD_CONSTEXPR qd_real &operator*=(double a);
  QD_CONSTEXPR qd_real &operator*=(const dd_real &a);
  QD_CONSTEXPR qd_real &operator*=(const qd_real &a);

  QD_CONSTEXPR static qd_real sloppy_div(const qd_real &a, const dd_real &b);
  QD_CONSTEXPR static qd_real accurate_div(const qd_real &a, const dd_real &b);
  QD_CONSTEXPR static qd_real sloppy_div(const qd_real &a, const qd_real &b);
  QD_CONSTEXPR static qd_real accurate_div(const qd_real &a, const qd_real &b);

  constexpr qd_real &operator/=(double a);
  constexpr qd_real &operator/=(const dd_real &a);
  constexpr qd_real &operator/=(const qd_real &a);

  qd_real operator^(int n) const;

  constexpr qd_real operator+() const;
  constexpr qd_real operator-() const;

  constexpr qd_real &operator=(double a);
  constexpr qd_real &operator=(const dd_real &a);
#if 0
  qd_real &operator=(const char *s);
#endif

  QD_CONSTEXPR bool is_zero() const;
  QD_CONSTEXPR bool is_one() const;
  QD_CONSTEXPR bool is_positive() const;
  QD_CONSTEXPR bool is_negative() const;

  static qd_real rand(void);

  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int len, int precision = _ndigits, 
      bool showpos = false, bool uppercase = false) const;
  std::string to_strig(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;
  static int read(const char *s, qd_real &a);

  /* Debugging methods */
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "", 
                 std::ostream &os = std::cerr) const;

  static qd_real debug_rand();
  std::string to_string(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;

  //HL:
  operator dd_real() {
	  if(this->isnan()) return dd_real::_nan;
	  return dd_real(x[0],x[1]);
  }
};

namespace std {
  template <>
  class numeric_limits<qd_real> : public numeric_limits<double> {
  public:
    static QD_CONSTEXPR double epsilon() { return qd_real::_eps; }
    static QD_CONSTEXPR double (min)() { return qd_real::_min_normalized; }
    static QD_CONSTEXPR qd_real (max)() { return qd_real::_max; }
    static QD_CONSTEXPR qd_real safe_max() { return qd_real::_safe_max; }
    static const int digits = 209;
    static const int digits10 = 62;
  };
}


QD_API qd_real polyeval(const qd_real *c, int n, const qd_real &x);
QD_API qd_real polyroot(const qd_real *c, int n, 
    const qd_real &x0, int max_iter = 64, double thresh = 0.0);

QD_API qd_real qdrand(void);
QD_CONSTEXPR qd_real sqrt(const qd_real &a);

/* Computes  qd * d  where d is known to be a power of 2.
   This can be done component wise.                      */
QD_API QD_CONSTEXPR qd_real mul_pwr2(const qd_real &qd, double d);

QD_API QD_CONSTEXPR qd_real operator+(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR qd_real operator+(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR qd_real operator+(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR qd_real operator+(const qd_real &a, double b);
QD_API QD_CONSTEXPR qd_real operator+(double a, const qd_real &b);

QD_API QD_CONSTEXPR qd_real operator-(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR qd_real operator-(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR qd_real operator-(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR qd_real operator-(const qd_real &a, double b);
QD_API QD_CONSTEXPR qd_real operator-(double a, const qd_real &b);

QD_API QD_CONSTEXPR qd_real operator*(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR qd_real operator*(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR qd_real operator*(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR qd_real operator*(const qd_real &a, double b);
QD_API QD_CONSTEXPR qd_real operator*(double a, const qd_real &b);

QD_API QD_CONSTEXPR qd_real operator/(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR qd_real operator/(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR qd_real operator/(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR qd_real operator/(const qd_real &a, double b);
QD_API QD_CONSTEXPR qd_real operator/(double a, const qd_real &b);

QD_API QD_CONSTEXPR qd_real sqr(const qd_real &a);
QD_CONSTEXPR qd_real sqrt(const qd_real &a);
QD_CONSTEXPR qd_real pow(const qd_real &a, int n);
QD_CONSTEXPR qd_real pow(const qd_real &a, const qd_real &b);
QD_CONSTEXPR qd_real npwr(const qd_real &a, int n);

qd_real nroot(const qd_real &a, int n);

QD_CONSTEXPR qd_real rem(const qd_real &a, const qd_real &b);
QD_CONSTEXPR qd_real drem(const qd_real &a, const qd_real &b);
QD_CONSTEXPR qd_real divrem(const qd_real &a, const qd_real &b, qd_real &r);

QD_CONSTEXPR dd_real to_dd_real(const qd_real &a);
QD_CONSTEXPR double  to_double(const qd_real &a);
QD_CONSTEXPR int     to_int(const qd_real &a);

QD_API QD_CONSTEXPR bool operator==(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator==(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR bool operator==(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator==(double a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator==(const qd_real &a, double b);

QD_API QD_CONSTEXPR bool operator<(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator<(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR bool operator<(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator<(double a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator<(const qd_real &a, double b);

QD_API QD_CONSTEXPR bool operator>(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator>(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR bool operator>(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator>(double a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator>(const qd_real &a, double b);

QD_API QD_CONSTEXPR bool operator<=(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator<=(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR bool operator<=(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator<=(double a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator<=(const qd_real &a, double b);

QD_API QD_CONSTEXPR bool operator>=(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator>=(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR bool operator>=(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator>=(double a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator>=(const qd_real &a, double b);

QD_API QD_CONSTEXPR bool operator!=(const qd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator!=(const qd_real &a, const dd_real &b);
QD_API QD_CONSTEXPR bool operator!=(const dd_real &a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator!=(double a, const qd_real &b);
QD_API QD_CONSTEXPR bool operator!=(const qd_real &a, double b);

QD_CONSTEXPR qd_real fabs(const qd_real &a);
QD_CONSTEXPR qd_real abs(const qd_real &a);    /* same as fabs */

QD_CONSTEXPR qd_real ldexp(const qd_real &a, int n);

QD_CONSTEXPR qd_real nint(const qd_real& a);
QD_CONSTEXPR qd_real round(const qd_real& a);
QD_API QD_CONSTEXPR qd_real quick_nint(const qd_real &a);
QD_API qd_real floor(const qd_real &a);
QD_API qd_real ceil(const qd_real &a);
QD_API QD_CONSTEXPR qd_real aint(const qd_real &a);

QD_CONSTEXPR qd_real sin(const qd_real &a);
QD_CONSTEXPR qd_real cos(const qd_real &a);
QD_CONSTEXPR qd_real tan(const qd_real &a);
QD_CONSTEXPR void sincos(const qd_real &a, qd_real &s, qd_real &c);

QD_CONSTEXPR qd_real asin(const qd_real &a);
QD_CONSTEXPR qd_real acos(const qd_real &a);
QD_CONSTEXPR qd_real atan(const qd_real &a);
QD_CONSTEXPR qd_real atan2(const qd_real &y, const qd_real &x);

QD_CONSTEXPR qd_real exp(const qd_real &a);
QD_CONSTEXPR qd_real log(const qd_real &a);
QD_CONSTEXPR qd_real log10(const qd_real &a);

QD_CONSTEXPR qd_real sinh(const qd_real &a);
QD_CONSTEXPR qd_real cosh(const qd_real &a);
QD_CONSTEXPR qd_real tanh(const qd_real &a);
QD_CONSTEXPR void sincosh(const qd_real &a, qd_real &sin_qd, qd_real &cos_qd);

QD_CONSTEXPR qd_real asinh(const qd_real &a);
QD_CONSTEXPR qd_real acosh(const qd_real &a);
QD_CONSTEXPR qd_real atanh(const qd_real &a);

QD_API qd_real qdrand(void);

QD_CONSTEXPR qd_real (max)(const qd_real &a, const qd_real &b);
QD_CONSTEXPR qd_real (max)(const qd_real &a, const qd_real &b, const qd_real &c);
QD_CONSTEXPR qd_real (min)(const qd_real &a, const qd_real &b);
QD_CONSTEXPR qd_real (min)(const qd_real &a, const qd_real &b, const qd_real &c);

QD_CONSTEXPR qd_real fmod(const qd_real &a, const qd_real &b);

QD_API std::ostream &operator<<(std::ostream &s, const qd_real &a);
QD_API std::istream &operator>>(std::istream &s, qd_real &a);

/*  qd_inline.h  */
/*
 * include/qd_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains small functions (suitable for inlining) in the quad-double
 * arithmetic package.
 */
#ifndef _QD_QD_INLINE_H
#define _QD_QD_INLINE_H

#include <cmath>

/********** Constructors **********/
inline QD_CONSTEXPR qd_real::qd_real(double x0, double x1, double x2, double x3) : x() {
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
  x[3] = x3;
}

inline qd_real::qd_real(const double *xx) : x() {
  x[0] = xx[0];
  x[1] = xx[1];
  x[2] = xx[2];
  x[3] = xx[3];
}

inline QD_CONSTEXPR qd_real::qd_real(double x0) : x() {
  x[0] = x0;
  x[1] = x[2] = x[3] = 0.0;
}

inline QD_CONSTEXPR qd_real::qd_real() : x() {
	x[0] = 0.0; 
	x[1] = 0.0; 
	x[2] = 0.0; 
	x[3] = 0.0; 
}

inline QD_CONSTEXPR qd_real::qd_real(const dd_real &a) : x() {
  x[0] = a._hi();
  x[1] = a._lo();
  x[2] = x[3] = 0.0;
}

inline QD_CONSTEXPR qd_real::qd_real(int i) : x() {
  x[0] = static_cast<double>(i);
  x[1] = x[2] = x[3] = 0.0;
}

/********** Accessors **********/
inline QD_CONSTEXPR double qd_real::operator[](int i) const {
  return x[i];
}

inline QD_CONSTEXPR double &qd_real::operator[](int i) {
  return x[i];
}

inline QD_CONSTEXPR bool qd_real::isnan() const {
  return QD_ISNAN(x[0]) || QD_ISNAN(x[1]) || QD_ISNAN(x[2]) || QD_ISNAN(x[3]);
}

/********** Renormalization **********/
namespace qd {
inline QD_CONSTEXPR void quick_renorm(double &c0, double &c1,
                         double &c2, double &c3, double &c4) {
  double t0, t1, t2, t3;
  double s;
  s  = qd::quick_two_sum(c3, c4, t3);
  s  = qd::quick_two_sum(c2, s , t2);
  s  = qd::quick_two_sum(c1, s , t1);
  c0 = qd::quick_two_sum(c0, s , t0);

  s  = qd::quick_two_sum(t2, t3, t2);
  s  = qd::quick_two_sum(t1, s , t1);
  c1 = qd::quick_two_sum(t0, s , t0);

  s  = qd::quick_two_sum(t1, t2, t1);
  c2 = qd::quick_two_sum(t0, s , t0);
  
  c3 = t0 + t1;
}

inline QD_CONSTEXPR void renorm(double &c0, double &c1,
                   double &c2, double &c3) {
  double s0, s1, s2 = 0.0, s3 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = qd::quick_two_sum(c2, c3, c3);
  s0 = qd::quick_two_sum(c1, s0, c2);
  c0 = qd::quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;
  if (s1 != 0.0) {
    s1 = qd::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0)
      s2 = qd::quick_two_sum(s2, c3, s3);
    else
      s1 = qd::quick_two_sum(s1, c3, s2);
  } else {
    s0 = qd::quick_two_sum(s0, c2, s1);
    if (s1 != 0.0)
      s1 = qd::quick_two_sum(s1, c3, s2);
    else
      s0 = qd::quick_two_sum(s0, c3, s1);
  }

  c0 = s0;
  c1 = s1;
  c2 = s2;
  c3 = s3;
}

inline QD_CONSTEXPR void renorm(double &c0, double &c1,
                   double &c2, double &c3, double &c4) {
  double s0, s1, s2 = 0.0, s3 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = qd::quick_two_sum(c3, c4, c4);
  s0 = qd::quick_two_sum(c2, s0, c3);
  s0 = qd::quick_two_sum(c1, s0, c2);
  c0 = qd::quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;

  s0 = qd::quick_two_sum(c0, c1, s1);
  if (s1 != 0.0) {
    s1 = qd::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0) {
      s2 = qd::quick_two_sum(s2, c3, s3);
      if (s3 != 0.0)
        s3 += c4;
      else
        s2 += c4;
    } else {
      s1 = qd::quick_two_sum(s1, c3, s2);
      if (s2 != 0.0)
        s2 = qd::quick_two_sum(s2, c4, s3);
      else
        s1 = qd::quick_two_sum(s1, c4, s2);
    }
  } else {
    s0 = qd::quick_two_sum(s0, c2, s1);
    if (s1 != 0.0) {
      s1 = qd::quick_two_sum(s1, c3, s2);
      if (s2 != 0.0)
        s2 = qd::quick_two_sum(s2, c4, s3);
      else
        s1 = qd::quick_two_sum(s1, c4, s2);
    } else {
      s0 = qd::quick_two_sum(s0, c3, s1);
      if (s1 != 0.0)
        s1 = qd::quick_two_sum(s1, c4, s2);
      else
        s0 = qd::quick_two_sum(s0, c4, s1);
    }
  }

  c0 = s0;
  c1 = s1;
  c2 = s2;
  c3 = s3;
}
}

inline QD_CONSTEXPR  void qd_real::renorm() {
  qd::renorm(x[0], x[1], x[2], x[3]);
}

inline QD_CONSTEXPR void qd_real::renorm(double &e) {
  qd::renorm(x[0], x[1], x[2], x[3], e);
}


/********** Additions ************/
namespace qd {

inline QD_CONSTEXPR void three_sum(double &a, double &b, double &c) {
  double t1, t2, t3;
  t1 = qd::two_sum(a, b, t2);
  a  = qd::two_sum(c, t1, t3);
  b  = qd::two_sum(t2, t3, c);
}

inline QD_CONSTEXPR void three_sum2(double &a, double &b, double &c) {
  double t1, t2, t3;
  t1 = qd::two_sum(a, b, t2);
  a  = qd::two_sum(c, t1, t3);
  b = t2 + t3;
}

}

/* quad-double + double */
inline QD_CONSTEXPR qd_real operator+(const qd_real &a, double b) {
  double c0, c1, c2, c3;
  double e;

  c0 = qd::two_sum(a[0], b, e);
  c1 = qd::two_sum(a[1], e, e);
  c2 = qd::two_sum(a[2], e, e);
  c3 = qd::two_sum(a[3], e, e);

  qd::renorm(c0, c1, c2, c3, e);

  return qd_real(c0, c1, c2, c3);
}

/* quad-double + double-double */
inline QD_CONSTEXPR qd_real operator+(const qd_real &a, const dd_real &b) {

  double s0, s1, s2, s3;
  double t0, t1;

  s0 = qd::two_sum(a[0], b._hi(), t0);
  s1 = qd::two_sum(a[1], b._lo(), t1);

  s1 = qd::two_sum(s1, t0, t0);

  s2 = a[2];
  qd::three_sum(s2, t0, t1);

  s3 = qd::two_sum(t0, a[3], t0);
  t0 += t1;

  qd::renorm(s0, s1, s2, s3, t0);
  return qd_real(s0, s1, s2, s3);
}


/* double + quad-double */
inline QD_CONSTEXPR qd_real operator+(double a, const qd_real &b) {
  return (b + a);
}

/* double-double + quad-double */
inline QD_CONSTEXPR qd_real operator+(const dd_real &a, const qd_real &b) {
  return (b + a);
}

namespace qd {

/* s = quick_three_accum(a, b, c) adds c to the dd-pair (a, b).
 * If the result does not fit in two doubles, then the sum is 
 * output into s and (a,b) contains the remainder.  Otherwise
 * s is zero and (a,b) contains the sum. */
inline QD_CONSTEXPR double quick_three_accum(double &a, double &b, double c) {
  double s;
  bool za, zb;

  s = qd::two_sum(b, c, b);
  s = qd::two_sum(a, s, a);

  za = (a != 0.0);
  zb = (b != 0.0);

  if (za && zb)
    return s;

  if (!zb) {
    b = a;
    a = s;
  } else {
    a = s;
  }

  return 0.0;
}

}

inline QD_CONSTEXPR qd_real qd_real::ieee_add(const qd_real &a, const qd_real &b) {
  int i, j, k;
  double s, t;
  double u, v;   /* double-length accumulator */
  double x[4] = {0.0, 0.0, 0.0, 0.0};
  
  i = j = k = 0;
  if (fb::abs(a[i]) > fb::abs(b[j]))
    u = a[i++];
  else
    u = b[j++];
  if (fb::abs(a[i]) > fb::abs(b[j]))
    v = a[i++];
  else
    v = b[j++];

  u = qd::quick_two_sum(u, v, v);
  
  while (k < 4) {
    if (i >= 4 && j >= 4) {
      x[k] = u;
      if (k < 3)
        x[++k] = v;
      break;
    }

    if (i >= 4)
      t = b[j++];
    else if (j >= 4)
      t = a[i++];
    else if (fb::abs(a[i]) > fb::abs(b[j])) {
      t = a[i++];
    } else
      t = b[j++];

    s = qd::quick_three_accum(u, v, t);

    if (s != 0.0) {
      x[k++] = s;
    }
  }

  /* add the rest. */
  for (k = i; k < 4; k++)
    x[3] += a[k];
  for (k = j; k < 4; k++)
    x[3] += b[k];

  qd::renorm(x[0], x[1], x[2], x[3]);
  return qd_real(x[0], x[1], x[2], x[3]);
}

inline QD_CONSTEXPR qd_real qd_real::sloppy_add(const qd_real &a, const qd_real &b) {
  /*
  double s0, s1, s2, s3;
  double t0, t1, t2, t3;
  
  s0 = qd::two_sum(a[0], b[0], t0);
  s1 = qd::two_sum(a[1], b[1], t1);
  s2 = qd::two_sum(a[2], b[2], t2);
  s3 = qd::two_sum(a[3], b[3], t3);

  s1 = qd::two_sum(s1, t0, t0);
  qd::three_sum(s2, t0, t1);
  qd::three_sum2(s3, t0, t2);
  t0 = t0 + t1 + t3;

  qd::renorm(s0, s1, s2, s3, t0);
  return qd_real(s0, s1, s2, s3, t0);
  */

  /* Same as above, but addition re-organized to minimize
     data dependency ... unfortunately some compilers are
     not very smart to do this automatically */
  double s0, s1, s2, s3;
  double t0, t1, t2, t3;

  double v0, v1, v2, v3;
  double u0, u1, u2, u3;
  double w0, w1, w2, w3;

  s0 = a[0] + b[0];
  s1 = a[1] + b[1];
  s2 = a[2] + b[2];
  s3 = a[3] + b[3];

  v0 = s0 - a[0];
  v1 = s1 - a[1];
  v2 = s2 - a[2];
  v3 = s3 - a[3];

  u0 = s0 - v0;
  u1 = s1 - v1;
  u2 = s2 - v2;
  u3 = s3 - v3;

  w0 = a[0] - u0;
  w1 = a[1] - u1;
  w2 = a[2] - u2;
  w3 = a[3] - u3;

  u0 = b[0] - v0;
  u1 = b[1] - v1;
  u2 = b[2] - v2;
  u3 = b[3] - v3;

  t0 = w0 + u0;
  t1 = w1 + u1;
  t2 = w2 + u2;
  t3 = w3 + u3;

  s1 = qd::two_sum(s1, t0, t0);
  qd::three_sum(s2, t0, t1);
  qd::three_sum2(s3, t0, t2);
  t0 = t0 + t1 + t3;

  /* renormalize */
  qd::renorm(s0, s1, s2, s3, t0);
  return qd_real(s0, s1, s2, s3);
}

/* quad-double + quad-double */
inline QD_CONSTEXPR qd_real operator+(const qd_real &a, const qd_real &b) {
#ifndef QD_IEEE_ADD
  return qd_real::sloppy_add(a, b);
#else
  return qd_real::ieee_add(a, b);
#endif
}



/********** Self-Additions ************/
/* quad-double += double */
inline QD_CONSTEXPR qd_real &qd_real::operator+=(double a) {
  *this = *this + a;
  return *this;
}

/* quad-double += double-double */
inline QD_CONSTEXPR qd_real &qd_real::operator+=(const dd_real &a) {
  *this = *this + a;
  return *this;
}

/* quad-double += quad-double */
inline QD_CONSTEXPR qd_real &qd_real::operator+=(const qd_real &a) {
  *this = *this + a;
  return *this;
}

/********** Unary Minus **********/
inline QD_CONSTEXPR qd_real qd_real::operator-() const {
  return qd_real(-x[0], -x[1], -x[2], -x[3]);
}

/********** Subtractions **********/
inline QD_CONSTEXPR qd_real operator-(const qd_real &a, double b) {
  return (a + (-b));
}

inline QD_CONSTEXPR qd_real operator-(double a, const qd_real &b) {
  return (a + (-b));
}

inline QD_CONSTEXPR qd_real operator-(const qd_real &a, const dd_real &b) {
  return (a + (-b));
}

inline QD_CONSTEXPR qd_real operator-(const dd_real &a, const qd_real &b) {
  return (a + (-b));
}

inline QD_CONSTEXPR qd_real operator-(const qd_real &a, const qd_real &b) {
  return (a + (-b));
}

/********** Self-Subtractions **********/
inline QD_CONSTEXPR qd_real &qd_real::operator-=(double a) {
  return ((*this) += (-a));
}

inline QD_CONSTEXPR qd_real &qd_real::operator-=(const dd_real &a) {
  return ((*this) += (-a));
}

inline QD_CONSTEXPR qd_real &qd_real::operator-=(const qd_real &a) {
  return ((*this) += (-a));
}


inline QD_CONSTEXPR qd_real operator*(double a, const qd_real &b) {
  return (b * a);
}

inline QD_CONSTEXPR qd_real operator*(const dd_real &a, const qd_real &b) {
  return (b * a);
}

inline QD_CONSTEXPR qd_real mul_pwr2(const qd_real &a, double b) {
  return qd_real(a[0] * b, a[1] * b, a[2] * b, a[3] * b);
}

/********** Multiplications **********/
inline QD_CONSTEXPR qd_real operator*(const qd_real &a, double b) {
  double p0, p1, p2, p3;
  double q0, q1, q2;
  double s0, s1, s2, s3, s4;

  p0 = qd::two_prod(a[0], b, q0);
  p1 = qd::two_prod(a[1], b, q1);
  p2 = qd::two_prod(a[2], b, q2);
  p3 = a[3] * b;

  s0 = p0;

  s1 = qd::two_sum(q0, p1, s2);

  qd::three_sum(s2, q1, p2);

  qd::three_sum2(q1, q2, p3);
  s3 = q1;

  s4 = q2 + p2;

  qd::renorm(s0, s1, s2, s3, s4);
  return qd_real(s0, s1, s2, s3);

}

/* quad-double * double-double */
/* a0 * b0                        0
        a0 * b1                   1
        a1 * b0                   2
             a1 * b1              3
             a2 * b0              4
                  a2 * b1         5
                  a3 * b0         6
                       a3 * b1    7 */
inline QD_CONSTEXPR qd_real operator*(const qd_real &a, const dd_real &b) {
  double p0, p1, p2, p3, p4;
  double q0, q1, q2, q3, q4;
  double s0, s1, s2;
  double t0, t1;

  p0 = qd::two_prod(a[0], b._hi(), q0);
  p1 = qd::two_prod(a[0], b._lo(), q1);
  p2 = qd::two_prod(a[1], b._hi(), q2);
  p3 = qd::two_prod(a[1], b._lo(), q3);
  p4 = qd::two_prod(a[2], b._hi(), q4);
  
  qd::three_sum(p1, p2, q0);
  
  /* Five-Three-Sum */
  qd::three_sum(p2, p3, p4);
  q1 = qd::two_sum(q1, q2, q2);
  s0 = qd::two_sum(p2, q1, t0);
  s1 = qd::two_sum(p3, q2, t1);
  s1 = qd::two_sum(s1, t0, t0);
  s2 = t0 + t1 + p4;
  p2 = s0;

  p3 = a[2] * b._hi() + a[3] * b._lo() + q3 + q4;
  qd::three_sum2(p3, q0, s1);
  p4 = q0 + s2;

  qd::renorm(p0, p1, p2, p3, p4);
  return qd_real(p0, p1, p2, p3);
}

/* quad-double * quad-double */
/* a0 * b0                    0
        a0 * b1               1
        a1 * b0               2
             a0 * b2          3
             a1 * b1          4
             a2 * b0          5
                  a0 * b3     6
                  a1 * b2     7
                  a2 * b1     8
                  a3 * b0     9  */
inline QD_CONSTEXPR qd_real qd_real::sloppy_mul(const qd_real &a, const qd_real &b) {
  double p0, p1, p2, p3, p4, p5;
  double q0, q1, q2, q3, q4, q5;
  double t0, t1;
  double s0, s1, s2;

  p0 = qd::two_prod(a[0], b[0], q0);

  p1 = qd::two_prod(a[0], b[1], q1);
  p2 = qd::two_prod(a[1], b[0], q2);

  p3 = qd::two_prod(a[0], b[2], q3);
  p4 = qd::two_prod(a[1], b[1], q4);
  p5 = qd::two_prod(a[2], b[0], q5);

  /* Start Accumulation */
  qd::three_sum(p1, p2, q0);

  /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
  qd::three_sum(p2, q1, q2);
  qd::three_sum(p3, p4, p5);
  /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
  s0 = qd::two_sum(p2, p3, t0);
  s1 = qd::two_sum(q1, p4, t1);
  s2 = q2 + p5;
  s1 = qd::two_sum(s1, t0, t0);
  s2 += (t0 + t1);

  /* O(eps^3) order terms */
  s1 += a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0] + q0 + q3 + q4 + q5;
  qd::renorm(p0, p1, s0, s1, s2);
  return qd_real(p0, p1, s0, s1);
}

inline QD_CONSTEXPR qd_real qd_real::accurate_mul(const qd_real &a, const qd_real &b) {
  double p0, p1, p2, p3, p4, p5;
  double q0, q1, q2, q3, q4, q5;
  double p6, p7, p8, p9;
  double q6, q7, q8, q9;
  double r0, r1;
  double t0, t1;
  double s0, s1, s2;

  p0 = qd::two_prod(a[0], b[0], q0);

  p1 = qd::two_prod(a[0], b[1], q1);
  p2 = qd::two_prod(a[1], b[0], q2);

  p3 = qd::two_prod(a[0], b[2], q3);
  p4 = qd::two_prod(a[1], b[1], q4);
  p5 = qd::two_prod(a[2], b[0], q5);

  /* Start Accumulation */
  qd::three_sum(p1, p2, q0);

  /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
  qd::three_sum(p2, q1, q2);
  qd::three_sum(p3, p4, p5);
  /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
  s0 = qd::two_sum(p2, p3, t0);
  s1 = qd::two_sum(q1, p4, t1);
  s2 = q2 + p5;
  s1 = qd::two_sum(s1, t0, t0);
  s2 += (t0 + t1);

  /* O(eps^3) order terms */
  p6 = qd::two_prod(a[0], b[3], q6);
  p7 = qd::two_prod(a[1], b[2], q7);
  p8 = qd::two_prod(a[2], b[1], q8);
  p9 = qd::two_prod(a[3], b[0], q9);

  /* Nine-Two-Sum of q0, s1, q3, q4, q5, p6, p7, p8, p9. */
  q0 = qd::two_sum(q0, q3, q3);
  q4 = qd::two_sum(q4, q5, q5);
  p6 = qd::two_sum(p6, p7, p7);
  p8 = qd::two_sum(p8, p9, p9);
  /* Compute (t0, t1) = (q0, q3) + (q4, q5). */
  t0 = qd::two_sum(q0, q4, t1);
  t1 += (q3 + q5);
  /* Compute (r0, r1) = (p6, p7) + (p8, p9). */
  r0 = qd::two_sum(p6, p8, r1);
  r1 += (p7 + p9);
  /* Compute (q3, q4) = (t0, t1) + (r0, r1). */
  q3 = qd::two_sum(t0, r0, q4);
  q4 += (t1 + r1);
  /* Compute (t0, t1) = (q3, q4) + s1. */
  t0 = qd::two_sum(q3, s1, t1);
  t1 += q4;

  /* O(eps^4) terms -- Nine-One-Sum */
  t1 += a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + q6 + q7 + q8 + q9 + s2;

  qd::renorm(p0, p1, s0, t0, t1);
  return qd_real(p0, p1, s0, t0);
}

inline QD_CONSTEXPR qd_real operator*(const qd_real &a, const qd_real &b) {
#ifdef QD_SLOPPY_MUL
  return qd_real::sloppy_mul(a, b);
#else
  return qd_real::accurate_mul(a, b);
#endif
}

/* quad-double ^ 2  = (x0 + x1 + x2 + x3) ^ 2
                    = x0 ^ 2 + 2 x0 * x1 + (2 x0 * x2 + x1 ^ 2)
                               + (2 x0 * x3 + 2 x1 * x2)           */
inline QD_CONSTEXPR qd_real sqr(const qd_real &a) {
  double p0, p1, p2, p3, p4, p5;
  double q0, q1, q2, q3;
  double s0, s1;
  double t0, t1;
  
  p0 = qd::two_sqr(a[0], q0);
  p1 = qd::two_prod(2.0 * a[0], a[1], q1);
  p2 = qd::two_prod(2.0 * a[0], a[2], q2);
  p3 = qd::two_sqr(a[1], q3);

  p1 = qd::two_sum(q0, p1, q0);

  q0 = qd::two_sum(q0, q1, q1);
  p2 = qd::two_sum(p2, p3, p3);

  s0 = qd::two_sum(q0, p2, t0);
  s1 = qd::two_sum(q1, p3, t1);

  s1 = qd::two_sum(s1, t0, t0);
  t0 += t1;

  s1 = qd::quick_two_sum(s1, t0, t0);
  p2 = qd::quick_two_sum(s0, s1, t1);
  p3 = qd::quick_two_sum(t1, t0, q0);

  p4 = 2.0 * a[0] * a[3];
  p5 = 2.0 * a[1] * a[2];

  p4 = qd::two_sum(p4, p5, p5);
  q2 = qd::two_sum(q2, q3, q3);

  t0 = qd::two_sum(p4, q2, t1);
  t1 = t1 + p5 + q3;

  p3 = qd::two_sum(p3, t0, p4);
  p4 = p4 + q0 + t1;

  qd::renorm(p0, p1, p2, p3, p4);
  return qd_real(p0, p1, p2, p3);

}

/********** Self-Multiplication **********/
/* quad-double *= double */
inline QD_CONSTEXPR qd_real &qd_real::operator*=(double a) {
  *this = (*this * a);
  return *this;
}

/* quad-double *= double-double */
inline QD_CONSTEXPR qd_real &qd_real::operator*=(const dd_real &a) {
  *this = (*this * a);
  return *this;
}

/* quad-double *= quad-double */
inline QD_CONSTEXPR qd_real &qd_real::operator*=(const qd_real &a) {
  *this = *this * a;
  return *this;
}

inline QD_CONSTEXPR qd_real operator/ (const qd_real &a, const dd_real &b) {
#ifdef QD_SLOPPY_DIV
  return qd_real::sloppy_div(a, b);
#else
  return qd_real::accurate_div(a, b);
#endif
}

inline QD_CONSTEXPR qd_real operator/(const qd_real &a, const qd_real &b) {
#ifdef QD_SLOPPY_DIV
  return qd_real::sloppy_div(a, b);
#else
  return qd_real::accurate_div(a, b);
#endif
}

/* double / quad-double */
inline QD_CONSTEXPR qd_real operator/(double a, const qd_real &b) {
  return qd_real(a) / b;
}

/* double-double / quad-double */
inline QD_CONSTEXPR qd_real operator/(const dd_real &a, const qd_real &b) {
  return qd_real(a) / b;
}

/********** Self-Divisions **********/
/* quad-double /= double */
inline QD_CONSTEXPR qd_real &qd_real::operator/=(double a) {
  *this = (*this / a);
  return *this;
}

/* quad-double /= double-double */
inline QD_CONSTEXPR qd_real &qd_real::operator/=(const dd_real &a) {
  *this = (*this / a);
  return *this;
}

/* quad-double /= quad-double */
inline QD_CONSTEXPR qd_real &qd_real::operator/=(const qd_real &a) {
  *this = (*this / a);
  return *this;
}


/********** Exponentiation **********/
inline qd_real qd_real::operator^(int n) const {
  return pow(*this, n);
}

/********** Miscellaneous **********/
inline QD_CONSTEXPR qd_real abs(const qd_real &a) {
  return (a[0] < 0.0) ? -a : a;
}

inline QD_CONSTEXPR qd_real fabs(const qd_real &a) {
  return abs(a);
}

/* Quick version.  May be off by one when qd is very close
   to the middle of two integers.                         */
inline QD_CONSTEXPR qd_real quick_nint(const qd_real &a) {
    qd_real r = qd_real(fb::round(a[0]), fb::round(a[1]),
        fb::round(a[2]), fb::round(a[3]));
  r.renorm();
  return r;
}

/*********** Assignments ************/
/* quad-double = double */
inline QD_CONSTEXPR qd_real &qd_real::operator=(double a) {
  x[0] = a;
  x[1] = x[2] = x[3] = 0.0;
  return *this;
}

/* quad-double = double-double */
inline QD_CONSTEXPR qd_real &qd_real::operator=(const dd_real &a) {
  x[0] = a._hi();
  x[1] = a._lo();
  x[2] = x[3] = 0.0;
  return *this;
}

/********** Equality Comparison **********/
inline QD_CONSTEXPR bool operator==(const qd_real &a, double b) {
  return (a[0] == b && a[1] == 0.0 && a[2] == 0.0 && a[3] == 0.0);
}

inline QD_CONSTEXPR bool operator==(double a, const qd_real &b) {
  return (b == a);
}

inline QD_CONSTEXPR bool operator==(const qd_real &a, const dd_real &b) {
  return (a[0] == b._hi() && a[1] == b._lo() && 
          a[2] == 0.0 && a[3] == 0.0);
}

inline QD_CONSTEXPR bool operator==(const dd_real &a, const qd_real &b) {
  return (b == a);
}

inline QD_CONSTEXPR bool operator==(const qd_real &a, const qd_real &b) {
  return (a[0] == b[0] && a[1] == b[1] && 
          a[2] == b[2] && a[3] == b[3]);
}


/********** Less-Than Comparison ***********/
inline QD_CONSTEXPR bool operator<(const qd_real &a, double b) {
  return (a[0] < b || (a[0] == b && a[1] < 0.0));
}

inline QD_CONSTEXPR bool operator<(double a, const qd_real &b) {
  return (b > a);
}

inline QD_CONSTEXPR bool operator<(const qd_real &a, const dd_real &b) {
  return (a[0] < b._hi() || 
          (a[0] == b._hi() && (a[1] < b._lo() ||
                            (a[1] == b._lo() && a[2] < 0.0))));
}

inline QD_CONSTEXPR bool operator<(const dd_real &a, const qd_real &b) {
  return (b > a);
}

inline QD_CONSTEXPR bool operator<(const qd_real &a, const qd_real &b) {
  return (a[0] < b[0] ||
          (a[0] == b[0] && (a[1] < b[1] ||
                            (a[1] == b[1] && (a[2] < b[2] ||
                                              (a[2] == b[2] && a[3] < b[3]))))));
}

/********** Greater-Than Comparison ***********/
inline QD_CONSTEXPR bool operator>(const qd_real &a, double b) {
  return (a[0] > b || (a[0] == b && a[1] > 0.0));
}

inline QD_CONSTEXPR bool operator>(double a, const qd_real &b) {
  return (b < a);
}

inline QD_CONSTEXPR bool operator>(const qd_real &a, const dd_real &b) {
  return (a[0] > b._hi() || 
          (a[0] == b._hi() && (a[1] > b._lo() ||
                            (a[1] == b._lo() && a[2] > 0.0))));
}

inline QD_CONSTEXPR bool operator>(const dd_real &a, const qd_real &b) {
  return (b < a);
}

inline QD_CONSTEXPR bool operator>(const qd_real &a, const qd_real &b) {
  return (a[0] > b[0] ||
          (a[0] == b[0] && (a[1] > b[1] ||
                            (a[1] == b[1] && (a[2] > b[2] ||
                                              (a[2] == b[2] && a[3] > b[3]))))));
}


/********** Less-Than-Or-Equal-To Comparison **********/
inline QD_CONSTEXPR bool operator<=(const qd_real &a, double b) {
  return (a[0] < b || (a[0] == b && a[1] <= 0.0));
}

inline QD_CONSTEXPR bool operator<=(double a, const qd_real &b) {
  return (b >= a);
}

inline QD_CONSTEXPR bool operator<=(const qd_real &a, const dd_real &b) {
  return (a[0] < b._hi() || 
          (a[0] == b._hi() && (a[1] < b._lo() || 
                            (a[1] == b._lo() && a[2] <= 0.0))));
}

inline QD_CONSTEXPR bool operator<=(const dd_real &a, const qd_real &b) {
  return (b >= a);
}

inline QD_CONSTEXPR bool operator<=(const qd_real &a, const qd_real &b) {
  return (a[0] < b[0] || 
          (a[0] == b[0] && (a[1] < b[1] ||
                            (a[1] == b[1] && (a[2] < b[2] ||
                                              (a[2] == b[2] && a[3] <= b[3]))))));
}

/********** Greater-Than-Or-Equal-To Comparison **********/
inline QD_CONSTEXPR bool operator>=(const qd_real &a, double b) {
  return (a[0] > b || (a[0] == b && a[1] >= 0.0));
}

inline QD_CONSTEXPR bool operator>=(double a, const qd_real &b) {
  return (b <= a);
}

inline QD_CONSTEXPR bool operator>=(const qd_real &a, const dd_real &b) {
  return (a[0] > b._hi() || 
          (a[0] == b._hi() && (a[1] > b._lo() || 
                            (a[1] == b._lo() && a[2] >= 0.0))));
}

inline QD_CONSTEXPR bool operator>=(const dd_real &a, const qd_real &b) {
  return (b <= a);
}

inline QD_CONSTEXPR bool operator>=(const qd_real &a, const qd_real &b) {
  return (a[0] > b[0] || 
          (a[0] == b[0] && (a[1] > b[1] ||
                            (a[1] == b[1] && (a[2] > b[2] ||
                                              (a[2] == b[2] && a[3] >= b[3]))))));
}



/********** Not-Equal-To Comparison **********/
inline QD_CONSTEXPR bool operator!=(const qd_real &a, double b) {
  return !(a == b);
}

inline QD_CONSTEXPR bool operator!=(double a, const qd_real &b) {
  return !(a == b);
}

inline QD_CONSTEXPR bool operator!=(const qd_real &a, const dd_real &b) {
  return !(a == b);
}

inline QD_CONSTEXPR bool operator!=(const dd_real &a, const qd_real &b) {
  return !(a == b);
}

inline QD_CONSTEXPR bool operator!=(const qd_real &a, const qd_real &b) {
  return !(a == b);
}



inline QD_CONSTEXPR qd_real aint(const qd_real &a) {
  return (a[0] >= 0) ? floor(a) : ceil(a);
}

inline QD_CONSTEXPR bool qd_real::is_zero() const {
  return (x[0] == 0.0);
}

inline QD_CONSTEXPR bool qd_real::is_one() const {
  return (x[0] == 1.0 && x[1] == 0.0 && x[2] == 0.0 && x[3] == 0.0);
}

inline QD_CONSTEXPR bool qd_real::is_positive() const {
  return (x[0] > 0.0);
}

inline QD_CONSTEXPR bool qd_real::is_negative() const {
  return (x[0] < 0.0);
}

inline QD_CONSTEXPR dd_real to_dd_real(const qd_real &a) {
  return dd_real(a[0], a[1]);
}

inline QD_CONSTEXPR double to_double(const qd_real &a) {
  return a[0];
}

inline QD_CONSTEXPR int to_int(const qd_real &a) {
  return static_cast<int>(a[0]);
}

inline QD_CONSTEXPR qd_real inv(const qd_real &qd) {
  return 1.0 / qd;
}

inline QD_CONSTEXPR qd_real (max)(const qd_real &a, const qd_real &b) {
  return (a > b) ? a : b;
}

inline QD_CONSTEXPR qd_real (max)(const qd_real &a, const qd_real &b,
                   const qd_real &c) {
  return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

inline QD_CONSTEXPR qd_real (min)(const qd_real &a, const qd_real &b) {
  return (a < b) ? a : b;
}

inline QD_CONSTEXPR qd_real (min)(const qd_real &a, const qd_real &b,
                   const qd_real &c) {
  return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

/* Random number generator */
inline qd_real qd_real::rand() {
  return qdrand();
}

inline QD_CONSTEXPR qd_real ldexp(const qd_real &a, int n) {
  return qd_real(fb::ldexp(a[0], n), fb::ldexp(a[1], n), 
                 fb::ldexp(a[2], n), fb::ldexp(a[3], n));
}

#endif /* _QD_QD_INLINE_H */
/*  qd_const.inl.h  */
/*
 * src/qd_const.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Defines constants used in quad-double package.
 */

/* Some useful constants. */
inline QD_CONSTEXPR const qd_real qd_real::_2pi = qd_real(6.283185307179586232e+00,
                                      2.449293598294706414e-16,
                                      -5.989539619436679332e-33,
                                      2.224908441726730563e-49);
inline QD_CONSTEXPR const qd_real qd_real::_pi = qd_real(3.141592653589793116e+00,
                                     1.224646799147353207e-16,
                                     -2.994769809718339666e-33,
                                     1.112454220863365282e-49);
inline QD_CONSTEXPR const qd_real qd_real::_pi2 = qd_real(1.570796326794896558e+00,
                                      6.123233995736766036e-17,
                                      -1.497384904859169833e-33,
                                      5.562271104316826408e-50);
inline QD_CONSTEXPR const qd_real qd_real::_pi4 = qd_real(7.853981633974482790e-01,
                                      3.061616997868383018e-17,
                                      -7.486924524295849165e-34,
                                      2.781135552158413204e-50);
inline QD_CONSTEXPR const qd_real qd_real::_3pi4 = qd_real(2.356194490192344837e+00,
                                       9.1848509936051484375e-17,
                                       3.9168984647504003225e-33,
                                      -2.5867981632704860386e-49);
inline QD_CONSTEXPR const qd_real qd_real::_e = qd_real(2.718281828459045091e+00,
                                    1.445646891729250158e-16,
                                    -2.127717108038176765e-33,
                                    1.515630159841218954e-49);
inline QD_CONSTEXPR const qd_real qd_real::_log2 = qd_real(6.931471805599452862e-01,
                                       2.319046813846299558e-17,
                                       5.707708438416212066e-34,
                                       -3.582432210601811423e-50);
inline QD_CONSTEXPR const qd_real qd_real::_log10 = qd_real(2.302585092994045901e+00,
                                        -2.170756223382249351e-16,
                                        -9.984262454465776570e-33,
                                        -4.023357454450206379e-49);
inline QD_CONSTEXPR const qd_real qd_real::_nan = qd_real(qd::_d_nan, qd::_d_nan,
                                      qd::_d_nan, qd::_d_nan);
inline QD_CONSTEXPR const qd_real qd_real::_inf = qd_real(qd::_d_inf, qd::_d_inf,
                                      qd::_d_inf, qd::_d_inf);

inline QD_CONSTEXPR const double qd_real::_eps = 1.21543267145725e-63; // = 2^-209
inline QD_CONSTEXPR const double qd_real::_min_normalized = 1.6259745436952323e-260; // = 2^(-1022 + 3*53)
inline QD_CONSTEXPR const qd_real qd_real::_max = qd_real(
    1.79769313486231570815e+308, 9.97920154767359795037e+291, 
    5.53956966280111259858e+275, 3.07507889307840487279e+259);
inline QD_CONSTEXPR const qd_real qd_real::_safe_max = qd_real(
    1.7976931080746007281e+308,  9.97920154767359795037e+291, 
    5.53956966280111259858e+275, 3.07507889307840487279e+259);
inline QD_CONSTEXPR const int qd_real::_ndigits = 62;

#include <qd/qd_real.inl.h>

#endif /* _QD_QD_REAL_H */

