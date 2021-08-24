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
#include "qd_config.h"

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
#include "dd_const.inl.h"

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



#include "dd_inline.h"
#endif /* _QD_DD_REAL_H */

