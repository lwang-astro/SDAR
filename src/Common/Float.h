#pragma once

#include <limits>

#ifdef USE_QD
#include <qd/qd_real.h>
#include <qd/qd_inline.h>
typedef qd_real Float;
const Float ROUND_OFF_ERROR_LIMIT=1e-60;
const Float NUMERIC_FLOAT_MAX = std::numeric_limits<double>::max();;

#elif USE_DD
#include <qd/dd_real.h>
#include <qd/dd_inline.h>
typedef dd_real Float;
const Float ROUND_OFF_ERROR_LIMIT=1e-30;
const Float NUMERIC_FLOAT_MAX = std::numeric_limits<double>::max();

#else
#include <limits>
typedef double Float;
#define to_int(x)     int(x)
#define to_double(x)  double(x)
const Float ROUND_OFF_ERROR_LIMIT=1e-14;
const Float NUMERIC_FLOAT_MAX = std::numeric_limits<Float>::max();

using std::sqrt;
using std::abs;
using std::pow;
using std::atan2;
using std::sin;
using std::cos;
using std::isnan;
using std::isinf;
#endif

