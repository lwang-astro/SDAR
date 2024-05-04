#pragma once

#include <limits>

#ifdef USE_QD
#include <qd/qd_real.h>
#include <qd/qd_inline.h>
typedef qd_real Float;
const Float ROUND_OFF_ERROR_LIMIT=1e-60;
const Float NUMERIC_FLOAT_MAX = std::numeric_limits<double>::max();;
const int WRITE_WIDTH=68;
const int WRITE_PRECISION=60;

#elif USE_DD
#include <qd/dd_real.h>
#include <qd/dd_inline.h>
typedef dd_real Float;
const Float ROUND_OFF_ERROR_LIMIT=1e-30;
const Float NUMERIC_FLOAT_MAX = std::numeric_limits<double>::max();
const int WRITE_WIDTH=38;
const int WRITE_PRECISION=30;

#elif USE_MPFRC
#include <mpreal.h>
using mpfr::mpreal;
typedef mpreal Float;
#define to_int(x)     int(x)
const int DIGITS_MPFRC = 80;
const Float ROUND_OFF_ERROR_LIMIT = std::numeric_limits<mpreal>::epsilon();
const Float NUMERIC_FLOAT_MAX = std::numeric_limits<mpreal>::max();
const int WRITE_WIDTH = DIGITS_MPFRC + 8;
const int WRITE_PRECISION = DIGITS_MPFRC;

#else
#include <limits>
typedef double Float;
#define to_int(x)     int(x)
#define to_double(x)  double(x)
const Float ROUND_OFF_ERROR_LIMIT=1e-14;
const Float NUMERIC_FLOAT_MAX = std::numeric_limits<Float>::max();
const int WRITE_WIDTH=23;
const int WRITE_PRECISION=14;

using std::sqrt;
using std::abs;
using std::pow;
using std::atan2;
using std::acos;
using std::sin;
using std::cos;
//#ifndef __INTEL_COMPILER
//using std::isnan;
//using std::isinf;
//#endif
#endif

#if (defined USE_QD) || (defined USE_DD) || (defined USE_MPFRC)
#define ISNAN(x) isnan(x)
#define ISINF(x) isinf(x)
#else
#define ISNAN(x) std::isnan(x)
#define ISINF(x) std::isinf(x)
#endif
