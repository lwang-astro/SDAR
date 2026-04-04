#pragma once

#include <limits>
#include <cmath>

#ifdef USE_MPFRC
#include <mpreal.h>
using mpfr::mpreal;
typedef mpreal Float;
#define to_int(x)     int(x)
static Float ROUND_OFF_ERROR_LIMIT = std::numeric_limits<mpreal>::epsilon();
static Float NUMERIC_FLOAT_MAX = std::numeric_limits<mpreal>::max();
static Float FLOAT_NAN = std::numeric_limits<mpreal>::quiet_NaN();
static int WRITE_WIDTH = 38;
static int WRITE_PRECISION = 30;

inline long long int llroundFloat(const Float& x) {
    return static_cast<long long int>(floor(x + Float(0.5)));
}

//* set MPFRC digits
void setMPFRPrec(const int digits) {
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    ROUND_OFF_ERROR_LIMIT = std::numeric_limits<mpreal>::epsilon();
    NUMERIC_FLOAT_MAX = std::numeric_limits<mpreal>::max();
    FLOAT_NAN = std::numeric_limits<mpreal>::quiet_NaN();
    WRITE_WIDTH = digits + 8;
    WRITE_PRECISION = digits;
    std::cerr<<"Set MPFR precision to "<<digits<<" digits."<<std::endl;
}

#else
#include <limits>
typedef double Float;
#define to_int(x)     int(x)
#define to_double(x)  double(x)
const Float ROUND_OFF_ERROR_LIMIT=1e-14;
const Float NUMERIC_FLOAT_MAX = std::numeric_limits<Float>::max();
const Float FLOAT_NAN = std::numeric_limits<Float>::quiet_NaN();
const int WRITE_WIDTH=23;
const int WRITE_PRECISION=14;

inline long long int llroundFloat(const Float& x) {
    return std::llround(x);
}

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
