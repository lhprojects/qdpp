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
#include "qd_config.h"
#include "two_basics.h"
#include "double_basics.h"
#include "util.h"

QD_CONSTEXPR bool dd_real::isnan() const { return QD_ISNAN(x[0]) || QD_ISNAN(x[1]); }
QD_CONSTEXPR bool dd_real::isfinite() const { return QD_ISFINITE(x[0]); }
QD_CONSTEXPR bool dd_real::isinf() const { return QD_ISINF(x[0]); }

//#define QD_OP_TAIL(x) if(!std::is_constant_evaluated()) assert(!(x).isnan())
#define QD_OP_TAIL(x) do { } while(false)

/* Cast to double. */
inline QD_CONSTEXPR double to_double(const dd_real& a)
{
    return a.x[0];
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

inline QD_CONSTEXPR dd_real ldexp(const dd_real& a, int exp)
{
    return dd_real(fb::ldexp(a.x[0], exp), fb::ldexp(a.x[1], exp));
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

template<class It1, class It2, class Real>
inline QD_CONSTEXPR int qd_read(It1 first,
    It2 end,
    Real& md)
{
    int sign = 0;
    int point = -1;  /* location of decimal point */
    int nd = 0;      /* number of digits read */
    int32_t e = 0;       /* exponent. */
    int64_t e64 = 0;
    bool e_sign = false;
    bool done = false;
    bool has_int = false;
    bool has_frac = false;
    Real r = 0.0;  /* number being read */

    /* Skip any leading spaces */
    for (; !(done = (first == end)) && *first == ' '; ++first) {
    }

    // done set
    if (!done && *first == '-') {
        ++first;
        sign = -1;
    }

    for (; !(done = first == end); ++first) {
        char ch = *first;
        if (ch >= '0' && ch <= '9') {
            int d = ch - '0';
            r *= 10.0;
            r += (double)d;
            nd++;
            has_int = true;
        } else {
            break;
        }
    }
    if (!done && *first == '.') {
        ++first;
        point = nd;

        for (; !(done = first == end); ++first) {
            char ch = *first;
            if (ch >= '0' && ch <= '9') {
                // too many digits, will reduce accuracy
                if (nd <= sizeof(Real) / sizeof(double) * 19) {
                    int d = ch - '0';
                    r *= 10.0;
                    r += (double)d;
                    nd++;
                }
                has_frac = true;
            } else {
                break;
            }
        }
    }
    if (!has_int && !has_frac) {
        return -1;
    }

    if (!done && (*first == 'e' || *first == 'E')) {
        ++first;

        e64 = 0;
        e_sign = false;
        done = first == end;
        if (!done && *first == '+') {
            ++first;
        } else if (!done && *first == '-') {
            ++first;
            e_sign = true;
        }

        bool has_exp = false;
        for (; !(done = first == end); ++first) {
            char ch = *first;
            if (ch >= '0' && ch <= '9') {
                has_exp = true;
                int d = ch - '0';
                e64 *= 10;
                e64 += d;
                if (e64 > (int64_t(INT32_MAX) + e_sign)) {
                    // overflow
                    for (; !(first == end); ++first) {
                        char ch = *first;
                        if (!(ch >= '0' && ch <= '9')) {
                            break;
                        }
                    }
                    return -1;
                }
            } else {
                break;
            }
        }
        if (!has_exp) {
            return -1;
        }
    }


    if (e_sign) e = -e64;
    else e = e64;

    /* Adjust exponent to account for decimal point */
    if (point >= 0) {
        e -= (nd - point);
    }

    /* Multiply the the exponent */
    if (e != 0) {
        r *= (Real(10.0) ^ e);
    }

    md = (sign < 0) ? -r : r;
    return 0;
}

/* Reads in a double-double number from the string s. */
inline constexpr int dd_real::read(const char* s, dd_real& a)
{
    size_t len = std::char_traits<char>::length(s);
    return qd_read(s, s + len, a);
}

inline constexpr dd_real dd_real::read(const char* s)
{
    dd_real u;
    if (dd_real::read(s, u) < 0) {
        return dd_real::_nan;
    }
    return u;
}


namespace qd_literals {
    inline namespace dd {

        template<char... Chars>
        QD_CONSTEXPR dd_real to_dd_real()
        {
            char buffer[] = { Chars..., '\0' };
            return dd_real::read(buffer);
        }

        template<char... Chars>
        QD_CONSTEXPR dd_real operator ""_dd()
        {
            QD_CONSTEXPR dd_real u = to_dd_real<Chars...>();
            return u;
        }

#if defined(QD_USE_ULL_LITERAL)
        inline QD_CONSTEXPR dd_real operator""_dd(unsigned long long u)
        {
            if (u <= UINT64_MAX) {
                double hi = u >> 32;
                double lo = u & UINT32_MAX;
                hi = ::qd::two_sum(hi, lo, lo);
                return dd_real(hi, lo);
            } else {
                fprintf(stderr, "range overflow: can't convert unsigned long long"
                    "(%llu) to dd_real", u);
                std::terminate();
            }
        }
#endif
    }
}
#endif /* _QD_DD_INLINE_H */
