#pragma once

#ifdef USE_QD
#include <qd/dd_real.h>
#include <qd/dd_inline.h>
typedef dd_real Float;

#else
typedef double Float;
#define to_int(x)     int(x)
#define to_double(x)  double(x)

using std::sqrt;
using std::abs;
using std::pow;
using std::atan2;
using std::sin;
using std::cos;
#endif

