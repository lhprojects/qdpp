#ifndef QD_SPECIAL_FUNCTIONS_H
#define QD_SPECIAL_FUNCTIONS_H

#include <qd/qd_real.h>

qd_real tgamma(qd_real a);
qd_real tgamma(qd_real a, qd_real x);
qd_real tgamma_lower(qd_real a, qd_real x);

#include "special_functions.inl.h"

#endif
