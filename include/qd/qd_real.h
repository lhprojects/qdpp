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
#include "qd_config.h"
#include "dd_real.h"
#include "fpu.h"

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
  QD_CONSTEXPR qd_real(const dd_real& dd);
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

  QD_CONSTEXPR qd_real &operator/=(double a);
  QD_CONSTEXPR qd_real &operator/=(const dd_real &a);
  QD_CONSTEXPR qd_real &operator/=(const qd_real &a);

  QD_CONSTEXPR qd_real operator^(int n) const;

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

  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int len, int precision = _ndigits, 
      bool showpos = false, bool uppercase = false) const;
  std::string to_strig(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;
  QD_CONSTEXPR static int read(const char *s, qd_real &a);
  QD_CONSTEXPR static qd_real read(const char* s);

  /* Debugging methods */
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "", 
                 std::ostream &os = std::cerr) const;

  static qd_real debug_rand();
  std::string to_string(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;

  QD_CONSTEXPR explicit operator dd_real() {
      // we assume dd_real knows nothing about qd_real
      // so we use operator here
      if(this->isnan()) return dd_real::_nan;
	  return dd_real(x[0],x[1]);
  }
};


namespace qd_literals {
    inline namespace qd {
        QD_CONSTEXPR qd_real operator""_qd(char const* s);
#if defined(QD_USE_ULL_LITERAL)
        QD_CONSTEXPR qd_real operator""_qd(unsigned long long u);
#endif
    }
}

QD_API qd_real polyeval(const qd_real *c, int n, const qd_real &x);
QD_API qd_real polyroot(const qd_real *c, int n, 
    const qd_real &x0, int max_iter = 64, double thresh = 0.0);

// generate uniform random number in [0,1]
template<class Gen>
qd_real qdrand(Gen&);

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
QD_API QD_CONSTEXPR qd_real floor(const qd_real &a);
QD_API QD_CONSTEXPR qd_real ceil(const qd_real &a);
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

QD_CONSTEXPR qd_real (max)(const qd_real &a, const qd_real &b);
QD_CONSTEXPR qd_real (max)(const qd_real &a, const qd_real &b, const qd_real &c);
QD_CONSTEXPR qd_real (min)(const qd_real &a, const qd_real &b);
QD_CONSTEXPR qd_real (min)(const qd_real &a, const qd_real &b, const qd_real &c);

QD_CONSTEXPR qd_real fmod(const qd_real &a, const qd_real &b);

QD_API std::ostream &operator<<(std::ostream &s, const qd_real &a);
QD_API std::istream &operator>>(std::istream &s, qd_real &a);

#include "qd_inline.h"
#include "qd_const.inl.h"
#include "qd_real.inl.h"

// put limits at the tail to resolve dependecy
namespace std {
    template <>
    class numeric_limits<qd_real> : public numeric_limits<double> {
    public:
        static QD_CONSTEXPR double epsilon() { return qd_real::_eps; }
        static QD_CONSTEXPR double (min)() { return qd_real::_min_normalized; }
        static QD_CONSTEXPR qd_real(max)() { return qd_real::_max; }
        static QD_CONSTEXPR qd_real safe_max() { return qd_real::_safe_max; }
        static const int digits = 209;
        static const int digits10 = 62;
    };
}

#endif /* _QD_QD_REAL_H */

