#ifndef QD_FLOATING_MATH_H
#define QD_FLOATING_MATH_H

#include "double_basics.h"
#include "two_basics.h"
#include "dd_real.h"

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

        // calculate a * _2opi to [2 bits . 121 bits] is enough for double
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
            static_assert(1023 + 120 + B < (int)sizeof(_2opi));
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

    inline constexpr double exp(double a) noexcept
    {
        if (std::is_constant_evaluated()) {
            return exp_(a);
        } else {
            return std::exp(a);
        }
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
