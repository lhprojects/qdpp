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
    // TODO: currectly not accuracy
    inline constexpr double _ln_2max = 710.47586007394394204;
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

    inline constexpr double fabs_(double x)
    {
        if (std::is_constant_evaluated()) {
            return copysign(x, 0);
        } else {
            return std::fabs(x);
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
                    // dividend -> (divisor/2, divisor)
                    dividend_ -= divisor_;
                    // dividend -> (-divisor/2, 0)
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
