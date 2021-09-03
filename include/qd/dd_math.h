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

#include "qd_config.h"
#include "util.h"
#include "bits.h"
#include "double_math.h"
#include <vector>
#include <random>
#include <type_traits>

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

// a^b
// = exp(b log(a)) = exp(b log1p(a_))
inline QD_CONSTEXPR dd_real pow(const dd_real &a, const dd_real &b) {
  return exp(b * log(a));
}

inline QD_CONSTEXPR int n_inv_fact_dd = 17;
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
  { 2.81145725434552060e-15,  1.65088427308614326e-31},
  { 1.56192069685862252711E-16, 1.19106796602737540024E-32 },
  { 8.22063524662432949554E-18, 2.21418941196042653637E-34 },
};


inline QD_CONSTEXPR dd_real expn_tab[] = {
dd_real(2.71828182845904509080E+00,1.44564689172925015783E-16),
dd_real(7.38905609893065040694E+00,-1.79711394978391483333E-16),
dd_real(5.45981500331442362040E+01,2.87415780158441149725E-15),
dd_real(2.98095798704172830185E+03,-2.71032958168736330162E-14),
dd_real(8.88611052050787210464E+06,5.32118248350156400432E-10),
dd_real(7.89629601826806875000E+13,7.66097802263510790910E-03),
dd_real(6.23514908081161674391E+27,1.38997388724928466797E+11),
dd_real(3.88770840599459482143E+55,2.70796611036621685030E+39),
dd_real(1.51142766500410352771E+111,1.48059891676144569992E+94),
dd_real(2.28441358653975650489E+222,1.35492249440234442263E+206),
};
inline QD_CONSTEXPR dd_real expn_inv_tab[] = {
dd_real(3.67879441171442334024E-01,-1.24287536727883625847E-17),
dd_real(1.35335283236612702318E-01,-1.04238142328866898838E-17),
dd_real(1.83156388887341786686E-02,1.62506889942713985216E-18),
dd_real(3.35462627902511853224E-04,-1.44021825104257951948E-20),
dd_real(1.12535174719259116458E-07,-1.94396212385793011861E-24),
dd_real(1.26641655490941755372E-14,1.85890796267480905359E-31),
dd_real(1.60381089054863792659E-28,-7.36132522128421421118E-45),
dd_real(2.57220937264241481170E-56,1.51363675530534789420E-73),
dd_real(6.61626105670948527686E-112,-1.58342892691157677006E-129),
dd_real(4.37749103705305118351E-223,2.70607983022862282398E-239),
};
inline QD_CONSTEXPR dd_real exp_integer_dd(int a)
{
    dd_real f = 1.;
    if (a > 0) {
        for (int i = 0; i < 10 && a; ++i) {
            if (a & 1) {
                f *= expn_tab[i];
            }
            a >>= 1;
        }
    } else {
        a = -a;
        for (int i = 0; i < 10 && a; ++i) {
            if (a & 1) {
                f *= expn_inv_tab[i];
            }
            a >>= 1;
        }
    }
    return f;
}

/* Exponential.  Computes exp(x) in double-double precision. */
inline QD_CONSTEXPR dd_real exp_general(const dd_real &a, bool sub_one) {

  const double k = 512.0;
  const double inv_k = 1.0 / k;

  dd_real a_ = a;
  if (a.x[0] >= 709.0)
      return dd_real::_inf;

  if (a.x[0] <= -709.0)
    return 0.0;

  if (a.is_zero())
    return 1.0;

  if (a.is_one())
    return dd_real::_e;


  double a0 = a.x[0];
  int aint = (int)fb::round(a0);
  dd_real factor;
  if (aint) {
      a0 -= aint;
      a_.x[0] = qd::quick_two_sum(a0, a_.x[1], a_.x[1]);
      factor = exp_integer_dd(aint);
  }

  dd_real r = mul_pwr2(a_, inv_k);
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
  } while (fb::abs(to_double(t)) > inv_k * dd_real::_eps && i < 8);

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

  if (sub_one) {
      if (aint == 0) {
          return s;
      } else {
          return (1. + s) * factor - 1.;
      }
  } else {
      s += 1.0;
      if (aint == 0) {
          return s;
      } else {
          return s * factor;
      }
  }
}

inline QD_CONSTEXPR dd_real exp(const dd_real& a)
{
    return exp_general(a, false);
}
inline QD_CONSTEXPR dd_real expm1(const dd_real& a)
{
    return exp_general(a, true);
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
    if (fb::abs(a.x[0]) < 0.5) {
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
#if 1
  dd_real sq = sqrt((1. - a) * (1. + a));
  return atan2(a, sq);
#else
  return atan2(a, sqrt(1.0 - sqr(a)));
#endif
}

inline QD_CONSTEXPR dd_real acos(const dd_real &a) {
  dd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    return dd_real::_nan;
  }

  if (abs_a.is_one()) {
    return (a.is_positive()) ? dd_real(0.0) : dd_real::_pi;
  }
#if 1
  dd_real sq = sqrt((1. - a) * (1. + a));
  return atan2(sq, a);
#else
  return atan2(sqrt(1.0 - sqr(a)), a);
#endif
}
 
inline QD_CONSTEXPR dd_real sinh(const dd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  if (abs(a) > 0.5) {
    // exp(a) - exp(-a)
    dd_real ea = exp(a);
    return mul_pwr2(ea - inv(ea), 0.5);
  } else if (a.x[0] < -0.5) {
      dd_real ea = exp(-a);
      return mul_pwr2(inv(ea) - ea, 0.5);
  } else {
      dd_real ea = expm1(a);
      // 1 + ea - 1/(1+ea)
      // ea(2+ea)/(1+ea)
      return mul_pwr2(ea * (2. + ea) / (1. + ea), 0.5);
  }
}

inline dd_real QD_CONSTEXPR cosh(const dd_real &a) {
  if (a.is_zero()) {
    return 1.0;
  }
  if (a < 0.) {
      // when a << 0, ea is tiny and we may lose some precision
      dd_real ea = exp(-a);
      return mul_pwr2(ea + inv(ea), 0.5);
  } else {
      dd_real ea = exp(a);
      return mul_pwr2(ea + inv(ea), 0.5);
  }
}

dd_real QD_CONSTEXPR tanh(const dd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  if (a.x[0] >= 0.5) {
      // (1 - exp(-2a))/(1+exp(-2a))
      // 1-... = 2 exp(-2a)/(1+exp(-2a))
      dd_real ea = exp(-mul_pwr2(a, 2));
      return 1. - mul_pwr2(ea, 2) / (1 + ea);
  } else if (a.x[0] <= -0.5) {
      dd_real ea = exp(mul_pwr2(a, 2));
      return -1. + mul_pwr2(ea, 2) / (1 + ea);
  } else {
      // (exp(x) - exp(-x))/(exp(x) + exp(-x))
      // (ea + 1 - 1/(1 + ea))/(ea + 1 + 1/(1+ea))
      // ea*(2 + ea)/(1+sqr(1+ea))
      dd_real ea = expm1(a);
      return ea * (2. + ea) / (1. + sqr(1. + ea));
  }
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

inline QD_CONSTEXPR dd_real acosh(const dd_real &a) {
    if (a < 1.0) {
        return dd_real::_nan;
    }
#if 1
    if (a.x[0] < 1.25) {
        dd_real a_ = a - 1.;
        dd_real q = a_ * (a + 1.);
        dd_real sq = sqrt(q);
        return log1p(a_ + sq);
    } else {
        dd_real sq = sqrt((a - 1.) * (a + 1.));
        return log(a + sq);
    }
#elif 0
    dd_real sq = sqrt((a - 1.) * (a + 1.));
    return log(a + sq);
#elif 0
    return log(a + sqrt(sqr(a) - 1));
#endif
}

inline QD_CONSTEXPR dd_real atanh(const dd_real &a) {
    if (abs(a) >= 1.0) {
        return dd_real::_nan;
    }
    if (abs(a) < 0.25) {
        // (1.0 + a) / (1.0 - a) = 1 + r_
        dd_real r_ = mul_pwr2(a, 2.) / (1.0 - a);
        return mul_pwr2(log1p(r_), 0.5);
    } else {
        return mul_pwr2(log((1.0 + a) / (1.0 - a)), 0.5);
    }
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


inline constexpr dd_real ceil(const dd_real& a)
{
    double hi = fb::ceil(a.x[0]);
    double lo = 0.0;

    if (hi == a.x[0]) {
        /* High word is integer already.  Round the low word. */
        lo = fb::ceil(a.x[1]);
        hi = qd::quick_two_sum(hi, lo, lo);
    }

    return dd_real(hi, lo);
}

inline constexpr dd_real aint(const dd_real& a)
{
    return (a.x[0] >= 0.0) ? floor(a) : ceil(a);
}

template<class T>
inline QD_CONSTEXPR int qd_get_bits(T max_)
{
    int bits = 0;
    for (T m = max_; m;) {
        bits += 1;
        m >>= 1;
    }
    if (!((max_ + 1) & max_)) {
        bits += 1;
    }
    bits -= 1;
    return bits;
}

#ifdef QD_HAS_CONSTEXPR
static_assert(qd_get_bits(1) == 1);
static_assert(qd_get_bits(2) == 1);
static_assert(qd_get_bits(uint32_t(-1)) == 32);
static_assert(qd_get_bits(uint64_t(-1)) == 64);
static_assert(qd_get_bits(uint64_t(uint32_t(-1))) == 32);
static_assert(qd_get_bits(uint64_t(1) << 53) == 53);
#endif

/* generate uniform double in [0,1) */
template<class Gen>
inline double drand(Gen& gen)
{
    using result_type = typename Gen::result_type;
    const result_type max_ = (Gen::max)() - (Gen::min)();
    QD_CONSTEXPR int bits = qd_get_bits(max_);
    const result_type bits_mask = result_type(1) << (bits - 1);
    const int n_gen_p_d = (53 + bits - 1) / bits;
    const bool use_all_bits = !((max_ + 1) & max_);
    uint64_t u1 = 0;

    // compiler unroll the loop
    for (int i = 0; i < n_gen_p_d; ++i) {
        if (bits != 64)
            u1 <<= bits;
        result_type g;
        for (;;) {
            g = gen() - (Gen::min)();
            if (use_all_bits) {
                break;
            } else {
                if (g > bits_mask) {
                    if (i != 0) {
                        // not need first time
                        g &= bits_mask;
                    }
                    break;
                }
            }

        }
        u1 |= g;
    }

    // exact 53 bits? don't belive it, mask anyway
    const uint64_t df_mask = (uint64_t(1) << 53) - 1;
    u1 &= df_mask;

    double r1 = double(u1);
    r1 *= 1.110223024625156540E-16; // 2^{-53}
    return r1;


}
/* generate uniform double in [0,1), better quality */
template<class Gen>
inline double drand_fine(Gen& gen)
{
    const double base_f = 1.110223024625156540E-16;
    double r = drand(gen);
    if (r >= 0.5)
        return r;
    else {
        if (r != 0.) {
            double r2 = drand(gen) * base_f;
            // 0 <= r2 < r < 0.5
            double ra = r + r2;
            return ra;
        } else {
            double base = base_f;
            for (;;) {
                r = drand(gen);
                if (r != 0.) {
                    break;
                }
                base *= base_f;
            }
            r *= base;
            double r2 = drand(gen) * base * base_f;
            double ra = r + r2;
            return ra;
        }
    }
}

/* generate uniform dd_real in [0,1) */
template<class Gen>
inline dd_real ddrand(Gen &gen)
{
    double r1 = drand_fine(gen);
    double r2 = drand_fine(gen) * 1.110223024625156540E-16;
    return dd_real::add(r1, r2);
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

template<class T, class Gen>
std::enable_if_t<std::is_same_v<T, double>, double>
real_rand(Gen& gen)
{
    return drand(gen);
}

template<class T, class Gen>
std::enable_if_t<std::is_same_v<T, dd_real>, dd_real>
real_rand(Gen& gen)
{
    return ddrand(gen);
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

struct qd_std_rand_generator {
    using result_type = uint32_t;

    static QD_CONSTEXPR result_type max()
    {
        return RAND_MAX;
    }

    static QD_CONSTEXPR result_type min()
    {
        return 0;
    }

    result_type operator()() const
    {
        return std::rand();
    }
};

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




#endif

inline dd_real dd_real::debug_rand()
{

    if (std::rand() % 2 == 0) {
        qd_std_rand_generator gen;
        return ddrand(gen);
    }

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

