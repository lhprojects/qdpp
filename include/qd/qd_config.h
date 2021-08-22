/* include/qd/qd_config.h.  Generated from qd_config.h.in by configure.  */
#ifndef _QD_QD_CONFIG_H
#define _QD_QD_CONFIG_H  1

#include "config.h"

/* If fused multiply-add is available, define to correct macro for
   using it.  It is invoked as QD_FMA(a, b, c) to compute fl(a * b + c). 
   If correctly rounded multiply-add is not available (or if unsure), 
   keep it undefined.*/
#ifndef QD_FMA
/* #undef QD_FMA */
#endif

/* If fused multiply-subtract is available, define to correct macro for
   using it.  It is invoked as QD_FMS(a, b, c) to compute fl(a * b - c). 
   If correctly rounded multiply-subtract is not available (or if unsure), 
   keep it undefined.*/
#ifndef QD_FMS
/* #undef QD_FMS */
#endif

/* Set the following to 1 to make the addition and subtraction
   to satisfy the IEEE-style error bound 

      fl(a + b) = (1 + d) * (a + b)

   where |d| <= eps.  If set to 0, the addition and subtraction
   will satisfy the weaker Cray-style error bound

      fl(a + b) = (1 + d1) * a + (1 + d2) * b

   where |d1| <= eps and |d2| eps.           */
#ifndef QD_IEEE_ADD
#define QD_IEEE_ADD 1
#endif

/* Set the following to 1 to use slightly inaccurate but faster
   version of multiplication. */
// #define QD_SLOPPY_MUL 1

/* Set the following to 1 to use slightly inaccurate but faster
   version of division. */
// #define QD_SLOPPY_DIV 1

/* Define this macro to be the isfinite(x) function. */
#ifndef QD_ISFINITE
#define QD_ISFINITE(x) fb::isfinite(x)
#endif

/* Define this macro to be the isinf(x) function. */
#ifndef QD_ISINF
#define QD_ISINF(x) fb::isinf(x)
#endif

/* Define this macro to be the isnan(x) function. */
#ifndef QD_ISNAN
#define QD_ISNAN(x) fb::isnan(x)

#endif

/* For C++ only, inline constexpr */
#define QD_CONSTEXPR constexpr
#define QD_HAS_CONSTEXPR 1

#endif /* _QD_QD_CONFIG_H */
