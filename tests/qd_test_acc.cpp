#include <qd/dd.h>
#include <qd/qd_real.h>

#include "mpreal.h"
#include "mpreal2.h"

mpfr::mpreal get_mpreal(mpfr2::mpreal const &r)
{
    mpfr::mpreal a(r.mpfr_ptr());
    a.set_prec(mpfr::mpreal::get_default_prec(), MPFR_RNDN);
	return a;
}

mpfr::mpreal get_mpreal(double r)
{
	mpfr::mpreal a(r);
	return a;
}

mpfr::mpreal get_mpreal(dd_real r)
{
	mpfr::mpreal a(r.x[0]);
	a += r.x[1];
	return a;
}

mpfr::mpreal get_mpreal(qd_real r)
{
	mpfr::mpreal a(r.x[0]);
	a += r.x[1];
	a += r.x[2];
	a += r.x[3];
	return a;
}

double to_double(mpfr2::mpreal const &r)
{
	return r.toDouble();
}

double get_ups(dd_real r)
{
	int a;
	(void)std::frexp(r.x[0], &a);
	// 0.1xx... 01yy... 2^a
	//    52x    52y
	// thus 1ups = 2^{-1-52-2-52-a} 2^a
	double v = ldexp(1., -1 - 52 - 2 - 52 + a);
	if (v == 0) return fb::_subnorm_min;
	return v;
}

double get_ups(qd_real r)
{
	int a;
	(void)std::frexp(r.x[0], &a);
	// 0.1xx... 01yy... 01zz... 01ww... 2^a
	//    52x    52y    52z      52w
	// thus 1ups = 2^{-1-52-(2+52)*3-a} 2^a
	double v = ldexp(1., -1 - (52 + 2) * 3 - 52 + a);
	if (v == 0) return fb::_subnorm_min;
	return v;
}

double get_ups(double r)
{
	int a;
	(void)std::frexp(r, &a);
	// 0.1xx... 01yy... 01zz... 01ww... 2^a
	//    52x    52y    52z      52w
	// thus 1ups = 2^{-1-52-(2+52)*3-a} 2^a
	double v = ldexp(1., -1 - 52 + a);
	if (v == 0) return fb::_subnorm_min;
	return v;
}

double get_ups(mpfr2::mpreal const& r)
{
	double ups =  get_ups(r.toDouble(MPFR_RNDZ));
    ups = ldexp(ups, 53 - mpfr2::mpreal::get_default_prec());
	return ups;
}

template<class Real>
double get_err_ups(Real const &v, mpfr::mpreal const &ref) {
    mpfr::mpreal err = get_mpreal(v) - ref;
    double derr = err.toDouble();
    double ups = get_ups(v);
    double r = std::abs(derr / ups);
    return r;
}

gmp_randstate_t rstate;
int a = []() {
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate, 0); return 0;
}();

template<class T, class Gen>
std::enable_if_t<std::is_same_v<T, mpfr2::mpreal>, mpfr2::mpreal>
real_rand(Gen& gen)
{
	return mpfr2::urandom(rstate);
}

template<class Real, bool use_fb = false>
struct TestAcc {


	int64_t iters = 1000;
	double range = 10;
	template<class F>
	void test_accuracy(char const* com, F f)
	{
		std::mt19937_64 gen;

		int err_u[100] = {};
		double err_sum = 0;
		double err_cnt = 0;
		double err_max = 0;
		double a_max = 0;
		double b_max = 0;
		bool hasnan = false;
		std::uniform_int_distribution<int> uni(0, 1);
		for (int64_t i = 0; i < iters; ++i) {
			Real a = real_rand<Real>(gen) - 0.5;
			Real b = real_rand<Real>(gen) - 0.5;
			a = exp(range * a) * (uni(gen) ? 1 : -1);
			b = exp(range * b) * (uni(gen) ? 1 : -1);
			Real c = f(a, b);
			if (isnan(c)) {
				hasnan = true;
				f(a, b);
				continue;
			}
			mpfr::mpreal a_ = get_mpreal(a);
			mpfr::mpreal b_ = get_mpreal(b);
			mpfr::mpreal cref = f(a_, b_);
			double r = get_err_ups(c, cref);

			if (r < 0.5) {
				err_u[0] += 1;
			} else if (r < 1.) {
				err_u[1] += 1;
			} else if (r < 2.) {
				err_u[2] += 1;
			} else if (strcmp(com, "+") == 0 && r > 3) {
				err_u[3] += 1;
				Real c = f(a, b);
				Real d = f(a, b);
				Real e = f(a, b);
			}
			err_sum += r;
			err_cnt += 1;
			if (r > err_max) {
				a_max = to_double(a);
				b_max = to_double(b);
			}
			err_max = std::max(err_max, r);
		}
		printf("%15s %12.1f %12.1f %12.1f %12.1f %12.2f(%10.2E,%10.2E) %12s\n",
			com,
			100. * err_u[0] / err_cnt, 100. * err_u[1] / err_cnt,
			100. * err_u[2] / err_cnt,
			err_sum / err_cnt,
			err_max, a_max, b_max,
			hasnan ? "Y" : "N");
	}
	template<class F>
	void test_accuracy_uni(char const* com, F f)
	{
		std::mt19937_64 gen;

		int err_u[100] = {};
		double err_sum = 0;
		double err_cnt = 0;
		double err_max = 0;
		double a_max = 0;
		bool hasnan = false;
		std::uniform_int_distribution<int> uni(0, 1);
		for (int64_t i = 0; i < iters; ++i) {
			Real a = real_rand<Real>(gen) - 0.5;
			a = exp(range * a) * (uni(gen) ? 1 : -1);
			Real c = f(a);
			if (isnan(c)) {
				hasnan = true;
				f(a);
				continue;
			}
			mpfr::mpreal a_ = get_mpreal(a);
			mpfr::mpreal cref = f(a_);
			double r = get_err_ups(c, cref);
			if (r < 0.5) {
				err_u[0] += 1;
			} else if (r < 1.) {
				err_u[1] += 1;
			} else if (r < 2.) {
				err_u[2] += 1;
			} else if (com == std::string("sinh") && r > 100) {
				err_u[3] += 1;
				Real c = f(a);
				Real d = f(a);
				Real e = f(a);
			}
			err_sum += r;
			err_cnt += 1;
			if (r > err_max) {
				a_max = to_double(a);
			}
			err_max = std::max(err_max, r);
		}
		printf("%15s %12.1f %12.1f %12.1f %12.1f %12.2f(%21.2E) %12s\n",
			com,
			100. * err_u[0] / err_cnt, 100. * err_u[1] / err_cnt,
			100. * err_u[2] / err_cnt,
			err_sum / err_cnt,
			err_max, a_max,
			hasnan ? "Y" : "N");
	}

	TestAcc(double range, int64_t iters)
	{
		this->iters = iters;
		this->range = range;

		std::string name;
		if (std::is_same_v<Real,mpfr2::mpreal>) {
			name = "mpfr(prec=" + std::to_string(mpfr2::mpreal::get_default_prec()) +")";
		} else {
			name = sizeof(Real) == sizeof(dd_real) ?
				"dd_real" : (sizeof(Real) == sizeof(double) ? "double" : "qd_real");
		}
		printf("%s%s range=+-[%6.1E,%6.1E], iters=%7lld\n",
			name.c_str(),
			use_fb ? "(fb::*)" : "",
			exp(-0.5 * range), exp(0.5 * range),
			(long long)iters);
		printf("%15s %12s %12s %12s %12s %12s(%21s) %12s\n",
			"",
			"0-0.5ups[%]", "0.5-1ups[%]", "1-2ups[%]", "avg[ups]",
			"max[ups]", "x(,y)",
			"HasNan");

		{
			qd_real factor = exp_integer_qd(600);
			mpfr::mpreal ref = exp(mpfr::mpreal(600));
			double r = get_err_ups(factor, ref);
			assert(r < 1.);
		}

		test_accuracy("+", [](auto const &a, auto const& b) { return a + b; });
		test_accuracy("-", [](auto const& a, auto const& b) { return a - b; });
		test_accuracy("*", [](auto const& a, auto const& b) { return a * b; });
		test_accuracy("/", [](auto const &a, auto const& b) { return a / b; });

		if constexpr (std::is_same_v < Real, dd_real >
			|| std::is_same_v<Real, qd_real>) {
			test_accuracy("accurate_div", [](auto a, auto b) {
				using Type = std::remove_cv_t<decltype(a)>;
				if constexpr (std::is_same_v<Type, dd_real>
					|| std::is_same_v<Type, qd_real>) {
					return Type::accurate_div(a, b);
				} else {
					return a / b;
				}
				});
			test_accuracy("sloppy_div", [](auto a, auto b) {
				using Type = std::remove_cv_t<decltype(a)>;
				if constexpr (std::is_same_v<Type, dd_real>
					|| std::is_same_v<Type, qd_real>) {
					return Type::sloppy_div(a, b);
				} else {
					return a / b;
				}
				});

			test_accuracy("ieee_add", [](auto a, auto b) {
				using Type = std::remove_cv_t<decltype(a)>;
				if constexpr (std::is_same_v<Type, dd_real>
					|| std::is_same_v<Type, qd_real>) {
					return Type::ieee_add(a, b);
				} else {
					return a + b;
				}
				});
			test_accuracy("sloppy_add", [](auto a, auto b) {
				using Type = std::remove_cv_t<decltype(a)>;
				if constexpr (std::is_same_v<Type, dd_real>
					|| std::is_same_v<Type, qd_real>) {
					return Type::sloppy_add(a, b);
				} else {
					return a + b;
				}
				});
		}

		if (std::is_same_v<Real, qd_real>) {
			test_accuracy("ieee_add_fine", [](auto a, auto b) {
				using Type = std::remove_cv_t<decltype(a)>;
				if constexpr (std::is_same_v<Type, qd_real>) {
					return Type::ieee_add_fine(a, b);
				} else {
					return a + b;
				}
				});		
			test_accuracy("accurate_mul", [](auto a, auto b) {
				using Type = std::remove_cv_t<decltype(a)>;
				if constexpr (std::is_same_v<Type, qd_real>) {
					return Type::accurate_mul(a, b);
				} else {
					return a * b;
				}
				});
			test_accuracy("sloppy_mul", [](auto a, auto b) {
				using Type = std::remove_cv_t<decltype(a)>;
				if constexpr (std::is_same_v<Type, qd_real>) {
					return Type::sloppy_mul(a, b);
				} else {
					return a * b;
				}
				});
		}

		test_accuracy_uni("abs", [](auto const& a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::fabs_(a);
			}
			return fabs(a);
			});

		test_accuracy_uni("floor", [](auto a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::floor_(a);
			}
			return floor(a);
			});
		test_accuracy_uni("ceil", [](auto a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::ceil_(a);
			}
			return ceil(a);
			});
		test_accuracy_uni("round", [](auto a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::round_(a);
			}
			return round(a);
			});

		test_accuracy("fmod", [](auto a, auto b) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::fmod_(a, b);
			}
			return fmod(a, b);
			});


		test_accuracy_uni("sqrt", [](auto a) {
			if (a < 0.) return 0. * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::sqrt_(fb::fabs_(a));
			}
			return sqrt(a);
			});
		test_accuracy_uni("exp", [](auto a) {
			if (fabs(a) > 600.) return 0 * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::exp_(a);
			}
			return exp(a);
			});
		test_accuracy_uni("log", [](auto a) {
			if (a < 0.) return 0. * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::log_(fb::fabs_(a));
			}
			return log(a);
			});
		test_accuracy_uni("log10", [](auto a) {
			if (a < 0.) return 0. * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::log10_(fb::fabs_(a));
			}
			return log10(a);
			});
		test_accuracy("pow", [](auto a, auto b) {
			if (a < 0.) return 0. * a;
			if (abs(log(abs(a)) / fb::_d_ln2 * b) > 900) return 0. * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::pow_(fb::fabs_(a), b);
			}
			return pow(a, b);
			});


		test_accuracy_uni("sin", [](auto a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::sin_(a);
			}
			return sin(a);
			});
		test_accuracy_uni("cos", [](auto a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::cos_(a);
			}
			return cos(a);
			});
		test_accuracy_uni("tan", [](auto a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::tan_(a);
			}
			return tan(a);
			});
		test_accuracy_uni("asin", [](auto a) {
			if (a > 1. || a < -1.) return 0. * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::asin_(a);
			}
			return asin(a);
			});
		test_accuracy_uni("acos", [](auto a) {
			if (a > 1. || a < -1.) return 0. * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::acos_(a);
			}
			return acos(a);
			});
		test_accuracy_uni("atan", [](auto a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::atan_(a);
			}
			return atan(a);
			});

		test_accuracy_uni("sinh", [](auto a) {
			if (fabs(a) > 600.) return 0 * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::sinh_(a);
			}
			return sinh(a);
			});
		test_accuracy_uni("cosh", [](auto a) {
			if (fabs(a) > 600.) return 0 * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::cosh_(a);
			}
			return cosh(a);
			});
		test_accuracy_uni("tanh", [](auto a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::tanh_(a);
			}
			return tanh(a);
			});
		test_accuracy_uni("asinh", [](auto a) {
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::asinh_(a);
			}
			return asinh(a);
			});
		test_accuracy_uni("acosh", [](auto a) {
			if (a < 1.) return 0. * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::acosh_(a);
			}
			return acosh(a);
			});
		test_accuracy_uni("atanh", [](auto a) {
			if (a > 1. || a < -1.) return 0. * a;
			using Type = std::remove_cv_t<decltype(a)>;
			if constexpr (std::is_same_v<Type, double> && use_fb) {
				return fb::atanh_(a);
			}
			return atanh(a);
			});
	}

};
void test_accuracy_set()
{
	{
		double a = get_ups((dd_real(2.) / 3.).x[1]);
		double b = get_ups((dd_real(4.) / 5.).x[1]);
		double c = get_ups((dd_real(3.) / 5.).x[1]);
		double d = get_ups((dd_real(4.) / 7.).x[1]);
		double e = get_ups((dd_real(5.) / 7.).x[1]);
		double f = get_ups((dd_real(6.) / 7.).x[1]);
		double ups[] = { a,b,c,d,e,f };
		double g = *std::max_element(std::begin(ups), std::end(ups));
		assert(g == get_ups((dd_real(2.) / 3.)));
	}
	{
		double a = get_ups((qd_real(2.) / 3.).x[3]);
		double b = get_ups((qd_real(4.) / 5.).x[3]);
		double c = get_ups((qd_real(3.) / 5.).x[3]);
		double d = get_ups((qd_real(4.) / 7.).x[3]);
		double e = get_ups((qd_real(5.) / 7.).x[3]);
		double f = get_ups((qd_real(6.) / 7.).x[3]);
		double ups[] = { a,b,c,d,e,f };
		double g = *std::max_element(std::begin(ups), std::end(ups));
		assert(g == get_ups((qd_real(2.) / 3.)));

	}

	mpfr::mpreal::set_default_prec(500);
	int64_t n = 1'000;
#ifdef NDEBUG
	n = 100'000;
#endif


	int64_t f = 1;

	TestAcc<qd_real, false>(0.1, n * f);
	TestAcc<qd_real, false>(1, n * f);
	TestAcc<qd_real, false>(10, n * f);
	TestAcc<qd_real, false>(20, n * f);

	TestAcc<dd_real, false>(0.1, n * f);
	TestAcc<dd_real, false>(1, n * f);
	TestAcc<dd_real, false>(10, n * f);
	TestAcc<dd_real, false>(20, n * f);

	TestAcc<double, true>(30, n * f);
	TestAcc<double, true>(1, n * f);
	TestAcc<double, true>(10, n * f);
	TestAcc<double, true>(20, n * f);

	TestAcc<double, false>(0.1, n * f);
	TestAcc<double, false>(1, n * f);
	TestAcc<double, false>(10, n * f);
	TestAcc<double, false>(20, n * f);

	mpfr2::mpreal::set_default_prec(53);
	TestAcc<mpfr2::mpreal, false>(0.1, n * f);
	TestAcc<mpfr2::mpreal, false>(1, n * f);
	TestAcc<mpfr2::mpreal, false>(10, n * f);
	TestAcc<mpfr2::mpreal, false>(20, n * f);

	mpfr2::mpreal::set_default_prec(106);
	TestAcc<mpfr2::mpreal, false>(0.1, n * f);
	TestAcc<mpfr2::mpreal, false>(1, n * f);
	TestAcc<mpfr2::mpreal, false>(10, n * f);
	TestAcc<mpfr2::mpreal, false>(20, n * f);

	mpfr2::mpreal::set_default_prec(212);
	TestAcc<mpfr2::mpreal, false>(0.1, n * f);
	TestAcc<mpfr2::mpreal, false>(1, n * f);
	TestAcc<mpfr2::mpreal, false>(10, n * f);
	TestAcc<mpfr2::mpreal, false>(20, n * f);
}

int main()
{
	test_accuracy_set();
	return 0;
}
