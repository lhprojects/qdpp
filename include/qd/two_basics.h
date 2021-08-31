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
#include "qd_config.h"

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
#define QD_USE_MM_FMADD 1
#elif defined(__GNUC__) || defined(__clang__)
#if defined(__FMA__) || defined(__FMA4__)
#define QD_USE_MM_FMADD 1
#endif
#endif

#ifdef QD_USE_MM_FMADD
#include <immintrin.h>

    inline double QD_FMA_NOCHECK(double a, double b, double c)
    {
        __m128d aw = { a, 0 };
        __m128d bw = { b, 0 };
        __m128d cw = { c, 0 };
        __m128d answ = _mm_fmadd_sd(aw, bw, cw);
        return _mm_cvtsd_f64(answ);
    }

    inline double QD_FMS_NOCHECK(double a, double b, double c)
    {
        __m128d aw = { a, 0 };
        __m128d bw = { b, 0 };
        __m128d cw = { c, 0 };
        __m128d answ = _mm_fmsub_sd(aw, bw, cw);
        return _mm_cvtsd_f64(answ);
    }

#ifndef QD_FMA
#define QD_FMA QD_FMA_NOCHECK
#endif

#ifndef QD_FMS
#define QD_FMS QD_FMS_NOCHECK
#endif

#endif // #ifdef QD_USE_MM_FMADD

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
        if (!std::is_constant_evaluated()) {
            double p = a * b;
            err = QD_FMS(a, b, p);
            return p;
        } else
#endif
        {
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
        if (!std::is_constant_evaluated()) {
            double p = a * a;
            err = QD_FMS(a, a, p);
            return p;
        } else
#endif
        {
            double hi, lo;
            double q = a * a;
            split(a, hi, lo);
            err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
            return q;
        }

    }
}

#endif /* _QD_INLINE_H */
