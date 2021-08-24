/*
 * tests/qd_timer.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2004
 *
 * Contains code to time basic operations.
 */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <qd/qd_real.h>
#include <qd/fpu.h>
#include "tictoc.h"
#include <qd/dd.h>

using std::cout;
using std::cerr;
using std::endl;
using std::sqrt;
using std::strcmp;
using std::setw;
using std::setprecision;
using std::fixed;

// Global flags passed to the main program.
static bool flag_test_double = false;
static bool flag_test_dd = false;
static bool flag_test_qd = false;
static bool flag_verbose = false;
static int  long_factor = 1;

template <class T>
class TestSuite {
public:
  void test1();
  void test2();
  void test3();
  void test4();
  void test7();
  void testall();
  T pi();



  double factor() const
  {
      if (sizeof(T) == sizeof(double)) return 1.;
      else if (sizeof(T) == sizeof(dd_real)) return 1 / 30.;
      else if (sizeof(T) == sizeof(qd_real)) return 1 / 1000.;
      return 1;
  }
};

template <class T>
T TestSuite<T>::pi() { return T::_pi; }

template <>
double TestSuite<double>::pi() { return 3.141592653589793116; }

void print_timing(double nops, double t, char const *com) {
    double mops = 1.0e-6 * nops / t;
    printf("%25s:  %10.6f us  %10.4f MOP/s\n", com, 1./mops, mops);
}

template <class T>
void TestSuite<T>::test1() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing addition..." << endl;
  }

  //int n = 100000, i;
  int n = 100000000* factor(), i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 1.0 / T(7.0);
  T a2 = 1.0 / T(11.0);
  T a3 = 1.0 / T(13.0);
  T a4 = 1.0 / T(17.0);
  T b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    b1 += a1;
    b2 += a2;
    b3 += a3;
    b4 += a4;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << b1+b2+b3+b4 << endl;
    cout << 4*n << " operations in " << t << " s." << endl;
  }

  print_timing(4.0*n, t, "add");
}

template <class T>
void TestSuite<T>::test2() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing multiplication ..." << endl;
  }

  //int n = 100000, i;
  int n = 100000000 * factor(), i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 1.0 + 1.0 / T(static_cast<double>(n));
  T a2 = 1.0 + 2.0 / T(static_cast<double>(n));
  T a3 = 1.0 + 3.0 / T(static_cast<double>(n));
  T a4 = 1.0 + 4.0 / T(static_cast<double>(n));
  T b1 = 1.0, b2 = 1.0, b3 = 1.0, b4 = 1.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    b1 *= a1;
    b2 *= a2;
    b3 *= a3;
    b4 *= a4;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << b1+b2+b3+b4 << endl;
    cout << 4*n << " operations in " << t << " s." << endl;
  }

  print_timing(4.0*n, t, "mul");
}

template <class T>
void TestSuite<T>::test3() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing division ..." << endl;
  }

  //int n = 100000, i;
  int n = 10000000 * factor(), i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 1.0 + 1.0 / T(static_cast<double>(n));
  T a2 = 1.0 + 2.0 / T(static_cast<double>(n));
  T a3 = 1.0 + 3.0 / T(static_cast<double>(n));
  T a4 = 1.0 + 4.0 / T(static_cast<double>(n));
  T b1 = 1.0, b2 = 1.0, b3 = 1.0, b4 = 1.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    b1 /= a1;
    b2 /= a2;
    b3 /= a3;
    b4 /= a4;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << b1+b2+b3+b4 << endl;
    cout << 4*n << " operations in " << t << " s." << endl;
  }

  print_timing(4.0*n, t, "div");
}

template <class T>
void TestSuite<T>::test4() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing square root ..." << endl;
  }

  //int n = 10000, i;
  int n = 10000000 * factor(), i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
  T b1 = 1.0 + pi();
  T b2 = 2.0 + pi();
  T b3 = 3.0 + pi();
  T b4 = 4.0 + pi();

  tic(&tv);
  for (i = 0; i < n; i++) {
    a1 = sqrt(a1 + b1);
    a2 = sqrt(a2 + b2);
    a3 = sqrt(a3 + b3);
    a4 = sqrt(a4 + b4);
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << a1+a2+a3+a4 << endl;
    cout << 4*n << " operations in " << t << " s." << endl;
  }

  print_timing(4.0*n, t, "sqrt");
}


template <class T>
void TestSuite<T>::test7() {
  if (flag_verbose) {
    cout << endl;
    cout << "Timing dot ..." << endl;
  }

  //int n = 100000, i;
  int n = 1000000000 * factor(), i;
  tictoc tv;
  double t;
  n *= long_factor;

  T a1 = 1.0 / T(7.0);
  T a2 = 1.0 / T(11.0);
  T a3 = 1.0 / T(13.0);
  T a4 = 1.0 / T(17.0);
  T b1 = 1.0 - T(1.0) / static_cast<double>(n);
  T b2 = 1.0 - T(2.0) / static_cast<double>(n);
  T b3 = 1.0 - T(3.0) / static_cast<double>(n);
  T b4 = 1.0 - T(4.0) / static_cast<double>(n);
  T x1 = 1.0, x2 = 1.0, x3 = 1.0, x4 = 1.0;

  tic(&tv);
  for (i = 0; i < n; i++) {
    x1 = a1 + b1 * x1;
    x2 = a2 + b2 * x2;
    x3 = a3 + b3 * x3;
    x4 = a4 + b4 * x4;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << "n = " << n << "   t = " << t << endl;
    cout << "r = " << x1+x2+x3+x4 << endl;
    cout << 8*n << " operations in " << t << " s." << endl;
  }

  print_timing(8.0*n, t, "dot");
}
template<class F, class T>
void test_arg1(F f, T const &low, T const &hig, double n, char const * name)
{
    if (flag_verbose) {
        cout << endl;
        cout << "Timing " << name << " ..." << endl;
    }

    tictoc tv;

    T a = 0.0;
    T c = low;
    T d = (hig - low) / static_cast<double>(n);

    double t;
    tic(&tv);
    for (uint64_t i = 0; i < uint64_t(n); i++) {
        a = a + f(c);
        c += d;
    }

    t = toc(&tv);
    if (flag_verbose) {
        cout << "n = " << n << "   t = " << t << endl;
        cout << "a = " << a << endl;
        cout << n << " operations in " << t << " s." << endl;
    }

    print_timing(n, t, name);
}

template<class T, class F>
void test_arg0(F f, double n, char const* name)
{
    if (flag_verbose) {
        cout << endl;
        cout << "Timing " << name << " ..." << endl;
    }

    tictoc tv;

    T a = 0.0;

    double t;
    tic(&tv);
    for (uint64_t i = 0; i < uint64_t(n); i++) {
        a += f();
    }
    t = toc(&tv);

    if (flag_verbose) {
        cout << "n = " << n << "   t = " << t << endl;
        cout << "a = " << a << endl;
        cout << n << " operations in " << t << " s." << endl;
    }

    print_timing(n, t, name);
}

template <class T>
void TestSuite<T>::testall() {
  test1();
  test2();
  test3();
  test4();
  test7();

  double n = 1E8 * factor() * long_factor;

  test_arg1([](T const& a) ->T { return round(a);  }, T(-5.), T(5.), n, "round");
  test_arg1([](T const& a) ->T { return ceil(a);  }, T(-5.), T(5.), n, "ceil");
  test_arg1([](T const& a) ->T { return floor(a);  }, T(-5.), T(5.), n, "floor");


  test_arg1([](T const& a) ->T { return sin(a);  }, T(-5.), T(5.), n, "sin");
  test_arg1([](T const& a) ->T { return cos(a);  }, T(-5.), T(5.), n, "cos");
  test_arg1([](T const& a) ->T { return tan(a);  }, T(-5.), T(5.), n, "tan");
  
  test_arg1([](T const& a) ->T { return asin(a);  }, T(-1.), T(1.), n, "asin");
  test_arg1([](T const& a) ->T { return acos(a);  }, T(-1.), T(1.), n, "acos");
  test_arg1([](T const& a) ->T { return atan(a);  }, T(-5.), T(5.), n, "atan");
  
  test_arg1([](T const& a) ->T { return exp(a);  }, T(-5.), T(5.), n, "exp");

  test_arg1([](T const& a) ->T { return log(a);  }, T(0.1), T(5.), n, "log");
  test_arg1([](T const& a) ->T { return log10(a);  }, T(0.1), T(5.), n, "log10");

  test_arg1([](T const& a) ->T { return sinh(a);  }, T(-5.), T(5.), n, "sinh");
  test_arg1([](T const& a) ->T { return sinh(a);  }, T(-0.05), T(0.05), n, "sinh[-0.05,0.05]");
  test_arg1([](T const& a) ->T { return cosh(a);  }, T(-5.), T(5.), n, "cosh");
  test_arg1([](T const& a) ->T { return tanh(a);  }, T(-5.), T(5.), n, "tanh");
  test_arg1([](T const& a) ->T { return asinh(a);  }, T(-5.), T(5.), n, "asinh");
  test_arg1([](T const& a) ->T { return acosh(a);  }, T(1), T(5.), n, "acosh");
  test_arg1([](T const& a) ->T { return atanh(a);  }, T(-1), T(1.), n, "atanh");

  std::mt19937 gen32;
  std::mt19937_64 gen64;
  test_arg0<T>([&]() ->T { return real_rand<T>(gen32);  }, n, "real_rand[MT19937_32]");
  test_arg0<T>([&]() ->T { return real_rand<T>(gen64);  }, n, "real_rand[MT19937_64]");
  if (sizeof(T) == sizeof(double)) {
      test_arg0<T>([&]() ->double { return drand_fine(gen32);  }, n, "drand_fine[MT19937_32]");
      test_arg0<T>([&]() ->double { return drand_fine(gen64);  }, n, "drand_fine[MT19937_64]");
  }
}

void print_usage() {
  cout << "qd_test [-h] [-dd] [-qd] [-all]" << endl;
  cout << "  Performs timing tests of the quad-double library." << endl;
  cout << "  By default, double-double and quad-double arithmetics" << endl;
  cout << "  are timed." << endl;
  cout << endl;
  cout << "  -h -help  Prints this usage message." << endl;
  cout << "  -double   Time arithmetic of double." << endl;
  cout << "  -dd       Time arithmetic of double-double." << endl;
  cout << "  -qd       Time arithmetic of quad-double." << endl;
  cout << "  -all      Perform both double-double and quad-double tests." << endl;
  cout << "  -v        Verbose output." << endl;
  cout << "  -long     Perform a longer timing loop." << endl;
}

int main(int argc, char *argv[]) {
  unsigned int old_cw;
  fpu_fix_start(&old_cw);

#ifdef QD_FMA
  printf("FMA enabled\n");
#else
  printf("FMA disable\n");
#endif

  /* Parse the arguments. */
  char *arg;
  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "-help") == 0) {
      print_usage();
      std::exit(0);
    } else if (strcmp(arg, "-double") == 0) {
      flag_test_double = true;
    } else if (strcmp(arg, "-dd") == 0) {
      flag_test_dd = true;
    } else if (strcmp(arg, "-qd") == 0) {
      flag_test_qd = true;
    } else if (strcmp(arg, "-all") == 0) {
      flag_test_double = flag_test_dd = flag_test_qd = true;
    } else if (strcmp(arg, "-v") == 0) {
      flag_verbose = true;
    } else if (strcmp(arg, "-long") == 0) {
      long_factor *= 10;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }

  /* If no flag, test both double-double and quad-double. */
  if (!flag_test_double && !flag_test_dd && !flag_test_qd) {
	  flag_test_double=true;
    flag_test_dd = true;
    flag_test_qd = true;
  }

  if (flag_test_double) {
    TestSuite<double> test;

    cout << endl;
    cout << "Timing double" << endl;
    cout << "-------------" << endl;
    test.testall();
  }

  if (flag_test_dd) {
    TestSuite<dd_real> test;

    cout << endl;
    cout << "Timing dd_real" << endl;
    cout << "--------------" << endl;
    test.testall();
  }

  if (flag_test_qd) {
    TestSuite<qd_real> test;

    cout << endl;
    cout << "Timing qd_real" << endl;
    cout << "--------------" << endl;
    test.testall();
  }
  
  fpu_fix_end(&old_cw);
  return 0;
}

