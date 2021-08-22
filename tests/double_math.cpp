#include <qd/double_basics.h>
#include <qd/double_math.h>

int nerr = 0;
int neps = 2;
#define QdAssert(x) do {\
	bool const v = bool(x);\
	if (!v) {\
		++nerr;\
		char const *expr = #x;\
		printf("%2d: failed: %s\n", __LINE__, expr);\
	}\
} while(false)

double get_ups(double x)
{
	if (x == 0.) {
		return 0.;
	} else if(x > 0){
        if (std::nextafter(x, std::numeric_limits<double>::infinity())
			== std::numeric_limits<double>::infinity()) {
            return x - std::nextafter(x, 0);
        } else {
            return std::min(std::nextafter(x, std::numeric_limits<double>::infinity()) - x,
                x - std::nextafter(x, 0));
		}
	} else {
		return get_ups(-x);
	}
}
bool fb_close(double x, double y, int neps_ = neps)
{
	if (std::isnan(x) || std::isnan(y)) {
		return std::isnan(x) && std::isnan(y);
	} else if (std::isinf(x) || std::isinf(y)) {
		return x == y;
	} else {
		double ux = get_ups(x);
		double uy = get_ups(y);
		double u = std::max(ux, uy);
		double diff = x - y;
		bool v = std::abs(diff) <= neps_ * u;
		return v;
	}
}

#define QdClose(x, y) do {\
	bool v = fb_close(x, y);\
	if (!v) {\
		++nerr;\
		printf("%2d: failed: %s ~ %s\n", __LINE__, #x, #y);\
	}\
} while(false)

#define QdCloseEps(x, y, eps) do {\
	bool v = fb_close(x, y, eps);\
	if (!v) {\
		++nerr;\
		printf("%2d: failed: %s ~ %s (eps = %d)\n", __LINE__, #x, #y, (int)(eps));\
	}\
} while(false)

const int SUBNORM_MIN_EXP = -1074;
const int SUBNORM_MAX_EXP = -1023;
const int NORM_MIN_EXP = -1022;
const int NORM_MAX_EXP = 1023;

using namespace fb;

void test_isnan()
{
	static_assert(!fb::isinf_(std::numeric_limits<double>::quiet_NaN()));
	static_assert(!fb::isinf(0.0));
	static_assert(!fb::isinf(0.5));
	static_assert(fb::isinf(std::numeric_limits<double>::infinity()));
	static_assert(fb::isinf(-std::numeric_limits<double>::infinity()));

	static_assert(fb::isnan(std::numeric_limits<double>::quiet_NaN()));
	static_assert(!fb::isnan(0.0));
	static_assert(!fb::isnan(0.5));
	static_assert(!fb::isnan(std::numeric_limits<double>::infinity()));
	static_assert(!fb::isnan(-std::numeric_limits<double>::infinity()));

	static_assert(!fb::isfinite(std::numeric_limits<double>::quiet_NaN()));
	static_assert(fb::isfinite(0.0));
	static_assert(fb::isfinite(0.5));
	static_assert(!fb::isfinite(std::numeric_limits<double>::infinity()));
	static_assert(!fb::isfinite(-std::numeric_limits<double>::infinity()));
}

void test_signbit()
{
	static_assert(!fb::signbit(0.));
	static_assert(fb::signbit(-0.));
	QdAssert(!fb::signbit(0.));
	QdAssert(fb::signbit(-0.));
	static_assert(!fb::signbit(1.));
	static_assert(fb::signbit(-1.));
	QdAssert(!fb::signbit(1.));
	QdAssert(fb::signbit(-1.));
	static_assert(!fb::signbit(std::numeric_limits<double>::infinity()));
	static_assert(fb::signbit(-std::numeric_limits<double>::infinity()));
	QdAssert(!fb::signbit(std::numeric_limits<double>::infinity()));
	QdAssert(fb::signbit(-std::numeric_limits<double>::infinity()));

	static_assert(!fb::signbit(std::numeric_limits<double>::quiet_NaN()));
	static_assert(fb::signbit(-std::numeric_limits<double>::quiet_NaN()));
	QdAssert(!fb::signbit(std::numeric_limits<double>::quiet_NaN()));
	QdAssert(fb::signbit(-std::numeric_limits<double>::quiet_NaN()));
}

void test_abs()
{
	static_assert(std::bit_cast<uint64_t>(fb::abs(0.)) == std::bit_cast<uint64_t>(0.));
	static_assert(std::bit_cast<uint64_t>(fb::abs(-0.)) == std::bit_cast<uint64_t>(0.));
	QdAssert(std::bit_cast<uint64_t>(fb::abs(0.)) == std::bit_cast<uint64_t>(0.));
	QdAssert(std::bit_cast<uint64_t>(fb::abs(-0.)) == std::bit_cast<uint64_t>(0.));

	static_assert(std::bit_cast<uint64_t>(fb::abs(1.)) == std::bit_cast<uint64_t>(1.));
	static_assert(std::bit_cast<uint64_t>(fb::abs(-1.)) == std::bit_cast<uint64_t>(1.));
	QdAssert(std::bit_cast<uint64_t>(fb::abs(1.)) == std::bit_cast<uint64_t>(1.));
	QdAssert(std::bit_cast<uint64_t>(fb::abs(-1.)) == std::bit_cast<uint64_t>(1.));

	constexpr double inf = std::numeric_limits<double>::infinity();
	constexpr double nan = std::numeric_limits<double>::quiet_NaN();

	static_assert(std::bit_cast<uint64_t>(fb::abs(inf)) == std::bit_cast<uint64_t>(inf));
	static_assert(std::bit_cast<uint64_t>(fb::abs(-inf)) == std::bit_cast<uint64_t>(inf));
	QdAssert(std::bit_cast<uint64_t>(fb::abs(inf)) == std::bit_cast<uint64_t>(inf));
	QdAssert(std::bit_cast<uint64_t>(fb::abs(-inf)) == std::bit_cast<uint64_t>(inf));

	static_assert(std::bit_cast<uint64_t>(fb::abs(nan)) == std::bit_cast<uint64_t>(nan));
	static_assert(std::bit_cast<uint64_t>(fb::abs(-nan)) == std::bit_cast<uint64_t>(nan));
	QdAssert(std::bit_cast<uint64_t>(fb::abs(nan)) == std::bit_cast<uint64_t>(nan));
	QdAssert(std::bit_cast<uint64_t>(fb::abs(-nan)) == std::bit_cast<uint64_t>(nan));
}

void test_copysign()
{
	static_assert(std::bit_cast<uint64_t>(fb::copysign(0., 0.)) == std::bit_cast<uint64_t>(0.));
	static_assert(std::bit_cast<uint64_t>(fb::copysign(0., -0.)) == std::bit_cast<uint64_t>(-0.));
	QdAssert(std::bit_cast<uint64_t>(fb::copysign(0., 0.)) == std::bit_cast<uint64_t>(0.));
	QdAssert(std::bit_cast<uint64_t>(fb::copysign(0., -0.)) == std::bit_cast<uint64_t>(-0.));

}
void test_round()
{
	static_assert(std::bit_cast<uint64_t>(fb::round(0.)) == std::bit_cast<uint64_t>(0.));
	static_assert(std::bit_cast<uint64_t>(fb::round(-0.)) == std::bit_cast<uint64_t>(-0.));
	QdAssert(std::bit_cast<uint64_t>(fb::round(0.)) == std::bit_cast<uint64_t>(0.));
	QdAssert(std::bit_cast<uint64_t>(fb::round(-0.)) == std::bit_cast<uint64_t>(-0.));

	static_assert(fb::round(0.49999999999999994) == 0.);
	QdAssert(fb::round(0.49999999999999994) == 0.);
	static_assert(fb::round(-0.49999999999999994) == 0.);
	QdAssert(fb::round(-0.49999999999999994) == -0.);

	constexpr double inf = std::numeric_limits<double>::infinity();
	constexpr double nan = std::numeric_limits<double>::quiet_NaN();

	static_assert(fb::round(inf) == inf);
	QdAssert(fb::round(inf) == inf);
	static_assert(fb::round(-inf) == -inf);
	QdAssert(fb::round(-inf) == -inf);

	static_assert(fb::isnan(fb::round(nan)));
	QdAssert(fb::isnan(fb::round(-nan)));

	static_assert(fb::round(10.499999999999998) == 10.);
	QdAssert(fb::round(10.499999999999998) == 10.);
}

void test_floor()
{

	constexpr double inf = std::numeric_limits<double>::infinity();
	constexpr double nan = std::numeric_limits<double>::quiet_NaN();

	static_assert(fb::floor(inf) == inf);
	QdAssert(fb::floor(inf) == inf);
	static_assert(fb::floor(-inf) == -inf);
	QdAssert(fb::floor(-inf) == -inf);

	static_assert(fb::isnan(fb::floor(nan)));
	QdAssert(fb::isnan(fb::floor(-nan)));

	static_assert(std::bit_cast<uint64_t>(fb::floor(0.)) == std::bit_cast<uint64_t>(0.));
	static_assert(std::bit_cast<uint64_t>(fb::floor(-0.)) == std::bit_cast<uint64_t>(-0.));
	QdAssert(std::bit_cast<uint64_t>(fb::floor(0.)) == std::bit_cast<uint64_t>(0.));
	QdAssert(std::bit_cast<uint64_t>(fb::floor(-0.)) == std::bit_cast<uint64_t>(-0.));

	static_assert(fb::floor(0.5) == 0.0);
	QdAssert(fb::floor(0.5) == 0.0);

	static_assert(fb::floor(1.) == 1.);
	QdAssert(fb::floor(1.) == 1.);

	static_assert(fb::floor(1.5) == 1);
	QdAssert(fb::floor(1.5) == 1);

	static_assert(fb::floor(-0.5) == -1.0);
	QdAssert(fb::floor(-0.5) == -1.0);

	//static_assert(fb::floor(-1.) == -1.);
	QdAssert(fb::floor(-1.) == -1.);

	//static_assert(fb::floor(-1.5) == -2.);
	QdAssert(fb::floor(-1.5) == -2.);

	static_assert(fb::floor(_int_max) == _int_max);
	QdAssert(fb::floor(-_int_max) == -_int_max);

	static_assert(fb::floor(_max) == _max);
	QdAssert(fb::floor(-_max) == -_max);

}

void test_ceil()
{
	constexpr double inf = std::numeric_limits<double>::infinity();
	constexpr double nan = std::numeric_limits<double>::quiet_NaN();

	static_assert(fb::ceil(inf) == inf);
	QdAssert(fb::ceil(inf) == inf);
	static_assert(fb::ceil(-inf) == -inf);
	QdAssert(fb::ceil(-inf) == -inf);

	static_assert(fb::isnan(fb::ceil(nan)));
	QdAssert(fb::isnan(fb::ceil(-nan)));

	static_assert(std::bit_cast<uint64_t>(fb::ceil(0.)) == std::bit_cast<uint64_t>(0.));
	static_assert(std::bit_cast<uint64_t>(fb::ceil(-0.)) == std::bit_cast<uint64_t>(-0.));
	QdAssert(std::bit_cast<uint64_t>(fb::ceil(0.)) == std::bit_cast<uint64_t>(0.));
	QdAssert(std::bit_cast<uint64_t>(fb::ceil(-0.)) == std::bit_cast<uint64_t>(-0.));

	static_assert(fb::ceil(0.5) == 1.0);
	QdAssert(fb::ceil(0.5) == 1.0);

	static_assert(fb::ceil(1.) == 1.);
	QdAssert(fb::ceil(1.) == 1.);

	static_assert(fb::ceil(1.5) == 2.);
	QdAssert(fb::ceil(1.5) == 2.);

	static_assert(fb::ceil(-0.5) == -0.0);
	QdAssert(fb::ceil(-0.5) == -0.0);

	static_assert(fb::ceil(-1.) == -1.);
	QdAssert(fb::ceil(-1.) == -1.);

	static_assert(fb::ceil(-1.5) == -1.);
	QdAssert(fb::ceil(-1.5) == -1.);

	static_assert(fb::ceil(_int_max) == _int_max);
	QdAssert(fb::ceil(-_int_max) == -_int_max);

	static_assert(fb::ceil(_max) == _max);
	QdAssert(fb::ceil(-_max) == -_max);

}

void test_sqrt()
{
	constexpr double eps = std::numeric_limits<double>::epsilon();
	static_assert(fb::isnan(fb::sqrt_(-1.)));
	QdAssert(std::isnan(fb::sqrt_(-1.)));

	static_assert(fb::sqrt_(0.) == 0);
	QdAssert(fb::sqrt_(0.) == 0);

	QdAssert(std::abs(fb::sqrt_(_subnorm_min) / std::sqrt(_subnorm_min) - 1) <= 1 * eps);
	QdAssert(std::abs(fb::sqrt_(2.0 * _subnorm_min) / std::sqrt(2. * _subnorm_min) - 1) <= 1 * eps);
	QdAssert(std::abs(fb::sqrt_(9.0 * _subnorm_min) / std::sqrt(9. * _subnorm_min) - 1) <= 1 * eps);
	QdAssert(std::abs(fb::sqrt_(_subnorm_max) / std::sqrt(_subnorm_max) - 1) <= 1 * eps);


	QdAssert(fb::sqrt_(1 / 2.) == 0.7071067811865475244008443621048490392848359376884740365883398689);


	QdAssert(fb::sqrt_(1.) == 1.);
	QdAssert(fb::sqrt_(2.) == 1.4142135623730950488016887242096980785696718753769480731766797379);
	QdAssert(fb::sqrt_(3.) == 1.7320508075688772935274463415058723669428052538103806280558069794);
	QdAssert(fb::sqrt_(4.) == 2.);
	QdAssert(fb::sqrt_(5.) == 2.2360679774997896964091736687312762354406183596115257242708972454);
	QdAssert(fb::sqrt_(10.) == 3.1622776601683793319988935444327185337195551393252168268575048527);
	QdAssert(fb::sqrt_(100.) == 10.);
	QdAssert(fb::sqrt_(1000.) == 31.622776601683793319988935444327185337195551393252168268575048527);
	QdAssert(fb::sqrt_(10000.) == 100.0);
	QdAssert(fb::sqrt_(100000.) == 316.22776601683793319988935444327185337195551393252168268575048527);
	QdAssert(fb::sqrt_(1000000.) == 1000.);
	QdAssert(fb::sqrt_(10000000.) == 3162.2776601683793319988935444327185337195551393252168268575048527);
	QdAssert(fb::sqrt_(100000000.) == 10000.);
	QdAssert(fb::sqrt_(1000000000.) == 31622.776601683793319988935444327185337195551393252168268575048527);
	QdAssert(fb::sqrt_(1E10) == 1E5);
	QdAssert(fb::sqrt_(1E11) == 3.1622776601683793319988935444327185337195551393252168268575048527E5);
	QdAssert(fb::sqrt_(1E12) == 1E6);
	QdAssert(fb::sqrt_(1E13) == 3.1622776601683793319988935444327185337195551393252168268575048527E6);
	QdAssert(fb::sqrt_(1E14) == 1E7);
	QdAssert(fb::sqrt_(1E15) == 3.1622776601683793319988935444327185337195551393252168268575048527E7);

	// nextafeter(2, 0) = 2 - epsilon
	// (2 - 1/2 epsilon)^2 = 4 - 2 *epsilon + epsilon^2 > 4 - 2 *epsilon
	// so sqrt(4 - 2 *epsilon) =  2 - epsilon
	constexpr double m4 = 3.999999999999999555910790149937383830547332763671875;
	constexpr double sqrt_m4 = 1.9999999999999998889776975374843428761489221713404328034938459635;
	QdAssert(sqrt_m4 != 2.0);
	QdAssert(std::nextafter(4., 0) == 4 - 2 * eps);
	QdAssert(std::nextafter(2., 0) == 2 - eps);
	QdAssert(fb::sqrt_(std::nextafter(4., 0)) == sqrt_m4);
	QdAssert(std::sqrt(std::nextafter(4., 0)) == sqrt_m4);

	QdAssert(std::abs(fb::sqrt_(_int_max) / std::sqrt(_int_max) - 1) <= 1. * eps);
	QdAssert(std::abs(fb::sqrt_(_max) / std::sqrt(_max) - 1) <= 1. * eps);

}

void test_ldexp()
{
	static_assert(fb::ldexp_(0., 0) == 0.);
	QdAssert(fb::ldexp_(0., 0) == 0.);

	static_assert(fb::ldexp_(-0., 0) == 0.);
	QdAssert(fb::ldexp_(-0., 0) == 0.);

	static_assert(fb::ldexp_(-0., 1) == 0.);
	QdAssert(fb::ldexp_(-0., 1) == 0.);

	static_assert(fb::ldexp_(0.5, 1) == 1.);
	QdAssert(fb::ldexp_(0.5, 1) == 1.);

	static_assert(fb::ldexp_(1.0, 1) == 2.);
	QdAssert(fb::ldexp_(1.0, 1) == 2.);

	static_assert(fb::ldexp_(1.5, 1) == 3.);
	QdAssert(fb::ldexp_(1.5, 1) == 3.);

	
	double const test = 2 * _d_one_prev;
	QdAssert(std::ldexp(std::ldexp(test, EXP_NORMAL_MIN), -EXP_NORMAL_MIN) == test);
	QdAssert(fb::ldexp_(fb::ldexp_(test, EXP_NORMAL_MIN), -EXP_NORMAL_MIN) == test);

	// trunced
    QdAssert(std::ldexp(std::ldexp(test, EXP_NORMAL_MIN - 1), -EXP_NORMAL_MIN + 1) != test);
	QdAssert(fb::ldexp_(fb::ldexp_(test, EXP_NORMAL_MIN - 1), -EXP_NORMAL_MIN + 1) != test);

	QdAssert(fb::ldexp_(_subnorm_min, -SUBNORM_MIN_EXP) == 1.0);
	QdAssert(fb::ldexp_(_subnorm_min, -1) == 0);
	QdAssert(fb::ldexp_(_subnorm_min, 1) == 2 * _subnorm_min);

	//static_assert(fb::ldexp_(1.0, -1074) == _subnorm_min);
	QdAssert(fb::ldexp_(1.0, SUBNORM_MIN_EXP) == _subnorm_min);

	static_assert(fb::ldexp_(_d_one_prev, 1024) == _max);
	QdAssert(fb::ldexp_(_d_one_prev, 1024) == _max);

	static_assert(fb::ldexp_(_d_one, 1024) == _infinity);
	QdAssert(fb::ldexp_(_d_one, 1024) == _infinity);


	// round
	QdAssert(test * std::ldexp(1., EXP_NORMAL_MIN - 1) * std::ldexp(1., -(EXP_NORMAL_MIN - 1)) == 2.);
	QdAssert(fb::mul2pwr_(fb::mul2pwr_(test, EXP_NORMAL_MIN - 1), -(EXP_NORMAL_MIN - 1)) == 2.);

}
void test_frexp()
{
	auto frexp = [](double x) {
		int exp_;
		x = fb::frexp_(x, &exp_);
		return exp_;
	};
	QdAssert(frexp(1.) == 1);
	QdAssert(frexp(0.5) == 0);
	QdAssert(frexp(0.) == 0);
	QdAssert(frexp(-0.) == 0);
	QdAssert(frexp(_subnorm_min) == SUBNORM_MIN_EXP + 1);
	QdAssert(frexp(_subnorm_max) == SUBNORM_MAX_EXP + 1);
	QdAssert(frexp(_norm_min) == NORM_MIN_EXP + 1);
	QdAssert(frexp(_max) == NORM_MAX_EXP + 1);

}

void test_fmod()
{
	QdAssert(fb::fmod_(0., 2.) == 0.);
	QdAssert(fb::fmod_(-0., 2.) == 0.);
	QdAssert(fb::fmod_(1., 2.) == 1.);
	QdAssert(fb::fmod_(-1., 2.) == -1);
	QdAssert(fb::fmod_(0., -2.) == 0.);
	QdAssert(fb::fmod_(-0., -2.) == 0.);
	QdAssert(fb::fmod_(1., -2.) == 1.);
	QdAssert(fb::fmod_(-1., -2.) == -1.);

	QdAssert(fb::fmod_(0., 1.) == 0.);
	QdAssert(fb::fmod_(0., 2.) == 0.);
	QdAssert(fb::fmod_(1., 2.) == 1.);
	QdAssert(fb::fmod_(2., 2.) == 0.);
	QdAssert(fb::fmod_(3., 2.) == 1.);
	QdAssert(fb::fmod_(1., 3.) == 1.);
	QdAssert(fb::fmod_(2., 3.) == 2.);
	QdAssert(fb::fmod_(3., 3.) == 0.);
	QdAssert(fb::fmod_(4., 3.) == 1.);
	QdAssert(fb::fmod_(5., 3.) == 2.);
	QdAssert(fb::fmod_(6., 3.) == 0.);
	QdAssert(fb::fmod_(fb::ldexp(2, 100), 3.) == 2.);


	QdAssert(fb::fmod_(_subnorm_min, _subnorm_min) == 0);
	QdAssert(fb::fmod_(_subnorm_max, _subnorm_min) == 0);
	QdAssert(fb::fmod_(1., _subnorm_min) == 0);
	QdAssert(fb::fmod_(0., 2 * _subnorm_min) == 0.);
	QdAssert(fb::fmod_(_subnorm_min, 2 * _subnorm_min) == _subnorm_min);
	QdAssert(fb::fmod_(2 * _subnorm_min, 2 * _subnorm_min) == 0.);
	QdAssert(fb::fmod_(3 * _subnorm_min, 2 * _subnorm_min) == _subnorm_min);

}

void test_remainder()
{
	QdAssert(fb::remainder_(0., 2.) == 0.);
	QdAssert(fb::remainder_(-0., 2.) == 0.);
	QdAssert(fb::remainder_(1., 2.) == 1.);
	QdAssert(fb::remainder_(-1., 2.) == -1);
	QdAssert(fb::remainder_(0., -2.) == 0.);
	QdAssert(fb::remainder_(-0., -2.) == 0.);
	QdAssert(fb::remainder_(1., -2.) == 1.);
	QdAssert(fb::remainder_(-1., -2.) == -1.);

	QdAssert(fb::remainder_(0., 1.) == 0.);
	QdAssert(fb::remainder_(0., 2.) == 0.);

	QdAssert(fb::remainder_(1., 2.) == 1.);
	QdAssert(fb::remainder_(2., 2.) == 0.);
	QdAssert(fb::remainder_(3., 2.) == -1.);

	QdAssert(fb::remainder_(1., 3.) == 1.);
	QdAssert(fb::remainder_(2., 3.) == -1);
	QdAssert(fb::remainder_(3., 3.) == 0.);
	QdAssert(fb::remainder_(4., 3.) == 1.);
	QdAssert(fb::remainder_(5., 3.) == -1);
	QdAssert(fb::remainder_(6., 3.) == 0.);
	QdAssert(fb::remainder_(fb::ldexp(2, 100), 3.) == -1);


	QdAssert(fb::remainder_(_subnorm_min, _subnorm_min) == 0);
	QdAssert(std::remainder(_subnorm_min, _subnorm_min) == 0);
	QdAssert(fb::remainder_(_subnorm_max, _subnorm_min) == 0);
	QdAssert(std::remainder(_subnorm_max, _subnorm_min) == 0);
	QdAssert(fb::remainder_(1., _subnorm_min) == 0);
	QdAssert(std::remainder(1., _subnorm_min) == 0);
	QdAssert(fb::remainder_(0., 2 * _subnorm_min) == 0.);
	QdAssert(std::remainder(0., 2 * _subnorm_min) == 0.);
	QdAssert(fb::remainder_(_subnorm_min, 2 * _subnorm_min) == _subnorm_min);
	QdAssert(std::remainder(_subnorm_min, 2 * _subnorm_min) == _subnorm_min);
	QdAssert(fb::remainder_(2 * _subnorm_min, 2 * _subnorm_min) == 0.);
	QdAssert(std::remainder(2 * _subnorm_min, 2 * _subnorm_min) == 0.);
	QdAssert(fb::remainder_(3 * _subnorm_min, 2 * _subnorm_min) == -_subnorm_min);
	QdAssert(std::remainder(3 * _subnorm_min, 2 * _subnorm_min) == -_subnorm_min);

}

void test_log()
{
	QdAssert(std::isinf(fb::log_(0.)));
	QdAssert(std::isinf(fb::log_(1.)) == 0);

	neps = 1;
    QdClose(fb::log_(1. / 8), -2.079441541679835928251696364374529704226500403080765762362040028);
    QdClose(fb::log_(2. / 8), -1.386294361119890618834464242916353136151000268720510508241360018);
	QdClose(fb::log_(5. / 8), -0.470003629245735553650937031148342064700899048812248040449392137);
	QdClose(fb::log_(7. / 8), -0.133531392624522623146343620931349974589415673498904573902649878);
	QdClose(fb::log_(15. / 16), -0.064538521137571171672923915683992928128908625349753842835377812);

	QdAssert(std::abs(fb::log_(2.) / _d_ln2-1) <= fb::_eps);
	QdAssert(std::abs(fb::log_(3.) / _d_ln3 - 1) <= fb::_eps);
	QdAssert(std::abs(fb::log_(4.) / _d_ln4 - 1) <= fb::_eps);
	QdAssert(std::abs(fb::log_(5.) / _d_ln5 - 1) <= fb::_eps);
	QdAssert(std::abs(fb::log_(6.) / _d_ln6 - 1) <= fb::_eps);
	QdAssert(std::abs(fb::log_(7.) / _d_ln7 - 1) <= fb::_eps);
	QdAssert(std::abs(fb::log_(10.) / _d_ln10 - 1) <= fb::_eps);

	double ref;
	ref = 0.0307716586667536883710282075967721640916967399588903563498619953;
	QdAssert(std::abs(fb::log_(1 + std::ldexp(1, -5)) / ref - 1) <= fb::_eps);
	ref = 0.0009760859730554588959608249080171866726118343337845362377585982;
	QdAssert(std::abs(fb::log_(1 + std::ldexp(1, -10)) / ref - 1) <= fb::_eps);
	ref = 0.00000095367386165918823390841551496333614360314807097930285811225;
	QdAssert(std::abs(fb::log_(1 + std::ldexp(1, -20)) / ref - 1) <= fb::_eps);

	// log(pi trunced as double) = log(db::_pi);
	constexpr double _d_ln_d_pi = 1.144729885849400135161709026159303486216473844443;
	QdAssert(std::abs(fb::log_(fb::_pi) / _d_ln_d_pi - 1) <= fb::_eps);
}

void test_sin()
{
    QdAssert(fb::sin_(0.) == 0.);
    QdAssert(fb::sin_(_d_pi) == 1.2246467991473531772260659322749980E-16);
    QdAssert(fb::sin_(2 * _d_pi) == -2 * 1.2246467991473531772260659322749980);

	// test precision

    QdAssert(fb::sin_(1E-20) == 1E-20);
    QdClose(fb::sin_(1), 0.8414709848078965066525023216302989996225630607983710656727517099);
    QdClose(fb::sin_(1.5), 0.9974949866040544309417233711414873227066514259221158219499748240);
    QdClose(fb::sin_(0.875), 0.76754350223602703963457546705453980969304345480560);
	
	QdAssert(fb::sin_(_d_pi / 2) == 1.);
    QdAssert(fb::sin_(-_d_pi / 2) == -1.);

    QdAssert(fb::sin_(0.25 * _d_pi) > 0);
    QdAssert(fb::sin_(0.75 * _d_pi) > 0);
    QdAssert(fb::sin_(1.25 * _d_pi) < 0);
    QdAssert(fb::sin_(1.75 * _d_pi) < 0);

	// test huge
	QdClose(fb::sin_(1234567890), 0.9866539395014835605426723146285895945861167418886784440443422855);
	QdClose(std::sin(1234567890), 0.9866539395014835605426723146285895945861167418886784440443422855);
}

void test_cos()
{
	QdAssert(fb::cos_(0.) == 1.);
	QdAssert(fb::cos_(_d_pi) == -1.);
	QdAssert(fb::cos_(2 * _d_pi) == 1.);

	// test precision (0, pi/2)
	QdAssert(fb::cos_(1E-10) == 1.);
	QdClose(fb::cos_(0.875), 0.6409968581633251303565566227960341319230463941938503949089927381);
	QdClose(fb::cos_(1), 0.5403023058681397174009366074429766037323104206179222276700972553);
	QdClose(fb::cos_(1.5), 0.0707372016677029100881898514342687090850910275633468694226454171);
	QdClose(fb::cos_(1.5625), 0.0082962316238583774779064583372156289738830072478408249657);



	QdAssert(fb::cos_(_d_pi / 2) == 6.1232339957367658861303296613750015E-17);
	QdAssert(fb::cos_(-_d_pi / 2) == -6.1232339957367658861303296613750015E-17);

	QdAssert(fb::cos_(0.25 * _d_pi) > 0);
	QdAssert(fb::cos_(0.75 * _d_pi) < 0);
	QdAssert(fb::cos_(1.25 * _d_pi) < 0);
	QdAssert(fb::cos_(1.75 * _d_pi) > 0);

}



void test_exp()
{
	// 1 - ln(2) = 0.3068528194400546905827678785418234319244998656397447458793199905
	// printf("%.50f", 0.30685281944005471);
	// double(.) = 0.30685281944005471377323601700481958687305450439453
	constexpr double dx = -2.3190467971754404e-17;
	// exp(double(.))-1 = 0.35914091422952264919925780279754371593338985154563
	constexpr double a = 0.35914091422952264919925780279754371593338985154563;
    constexpr double e = 2 * (1. + a);
	double d = e - _d_e;
	QdClose(fb::exp_(0.), 1.);
    QdClose(fb::exp_(1.), _d_e);
	QdClose(fb::exp_(-1.), 1 / _d_e);
	QdClose(fb::exp_(2.), 7.3890560989306502272304274605750078131803155705518473240871278225);
	QdClose(fb::exp_(-2.), 0.1353352832366126918939994949724844034076315459095758814681588726);
}

void test_cbrt()
{
	QdAssert(fb::cbrt_(1) == 1.);
	QdAssert(fb::cbrt_(2) == 1.2599210498948731647672106072782283505702514647015079800819751121);
	QdAssert(fb::cbrt_(3) == 1.4422495703074083823216383107801095883918692534993505775464161945);
	QdAssert(fb::cbrt_(4) == 1.5874010519681994747517056392723082603914933278998530098082857618);
	QdAssert(fb::cbrt_(5) == 1.7099759466766969893531088725438601098680551105430549243828617074);
	QdAssert(fb::cbrt_(6) == 1.8171205928321396588912117563272605024282104631412196714813342979);
	QdAssert(fb::cbrt_(7) == 1.9129311827723891011991168395487602828624390503458757662106476404);
	QdAssert(fb::cbrt_(8) == 2.);
	QdAssert(fb::cbrt_(1000) == 10.);
	QdAssert(fb::cbrt_(0.001) == 0.1);
	QdAssert(fb::cbrt_(0.5625) == 0.8254818122236566709686524881022712392117211610796233342433);
	constexpr double b0 = fb::cbrt_(1);
	constexpr double b1 = fb::cbrt_(1000);
	constexpr double b2 = fb::cbrt_(2);
	constexpr double b3 = fb::cbrt_(0.001);
	constexpr double b4 = fb::cbrt_(0.5625);
}

void test_hypot()
{
	constexpr double a = fb::hypot_(0, 1);
	QdAssert(fb::hypot_(_infinity, 0) == _infinity);
	QdAssert(fb::hypot_(_infinity, 1) == _infinity);
	QdAssert(fb::hypot_(_infinity, _d_nan) == _infinity);
	QdAssert(fb::hypot_(-_infinity, 0) == _infinity);
	QdAssert(fb::hypot_(-_infinity, 1) == _infinity);
	QdAssert(fb::hypot_(-_infinity, _d_nan) == _infinity);

	QdAssert(fb::hypot_(0, 1) == 1);
    QdAssert(fb::hypot_(1, 0) == 1);
    QdAssert(fb::hypot_(1, 1) == 1.4142135623730950488016887242096980785696718753769480731766797379);
    QdAssert(fb::hypot_(1, 2) == 2.2360679774997896964091736687312762354406183596115257242708972454);
    QdAssert(fb::hypot_(1, 3) == 3.1622776601683793319988935444327185337195551393252168268575048527);
    QdAssert(fb::hypot_(1, 1000) == 1000.0004999998750000624999609375273437294922036132681579698943999);
    QdAssert(fb::hypot_(_d_pi, _d_pi / 2) == 3.51240736552036305965851275555178760679657218311525);
    QdAssert(fb::hypot_(_d_pi, _d_pi / 10) == 3.15726154208045479557128857463234333803672524107779);

	QdAssert(fb::hypot_(_infinity, 0, 1) == _infinity);
	QdAssert(fb::hypot_(_infinity, 1, 0) == _infinity);
	QdAssert(fb::hypot_(_infinity, _d_nan, 0) == _infinity);
	QdAssert(fb::hypot_(-_infinity, 0, 0) == _infinity);
	QdAssert(fb::hypot_(-_infinity, 1, 0) == _infinity);
	QdAssert(fb::hypot_(-_infinity, _d_nan, 0) == _infinity);

	QdAssert(fb::hypot_(0, 1, 0) == 1);
	QdAssert(fb::hypot_(1, 0, 1) == 1.4142135623730950488016887242096980785696718753769480731766797379);
	QdAssert(fb::hypot_(1, 0, 2) == 2.2360679774997896964091736687312762354406183596115257242708972454);
	QdAssert(fb::hypot_(1, 3, 0) == 3.1622776601683793319988935444327185337195551393252168268575048527);
	QdAssert(fb::hypot_(1, 1000, 0) == 1000.0004999998750000624999609375273437294922036132681579698943999);
	QdAssert(fb::hypot_(0, _d_pi, _d_pi / 2) == 3.51240736552036305965851275555178760679657218311525);
	QdAssert(fb::hypot_(0, _d_pi, _d_pi / 10) == 3.15726154208045479557128857463234333803672524107779);
}

void test_nan()
{
	constexpr double a = fb::nan_("1");
	double b = fb::nan_("1");
	double c = fb::nan_("");
}

int main()
{


	test_signbit();
	test_copysign();
	test_abs();
	test_ldexp();
	test_frexp();
	test_floor();
	test_round();
	test_ceil();
	test_sqrt();
	test_fmod();
	test_remainder();
	test_log();
	test_sin();
	test_cos();
	test_exp();
	test_cbrt();
	test_hypot();
	test_nan();

	if (nerr)
		return -1;
	else
		return 0;
}
