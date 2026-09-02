#pragma once
#include <boost/math/constants/constants.hpp>
#ifndef USE_MULTIPRECISION
using Real = double;
inline Real mp_sin(Real value) { return std::sin(value); }

inline Real mp_cos(Real value) { return std::cos(value); }

inline Real mp_exp(Real value) { return std::exp(value); }

inline Real mp_sqrt(Real value) { return std::sqrt(value); }

inline Real mp_pow(Real value, Real deg) { return std::pow(value, deg); }

inline Real mp_max(Real a, Real b) { return std::max(a, b); }
#endif

#ifdef USE_MULTIPRECISION

#include <boost/multiprecision/fwd.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost::multiprecision;
using Real = boost::multiprecision::static_mpfr_float_50;

inline Real mp_sin(Real value) { return boost::multiprecision::sin(value); }

inline Real mp_cos(Real value) { return boost::multiprecision::cos(value); }

inline Real mp_exp(Real value) { return exp(value); }

inline Real mp_sqrt(Real value) { return sqrt(value); }

inline Real mp_pow(Real value, Real deg) { return pow(value, deg); }
inline Real mp_max(Real a, Real b) { return a < b ? b : a; }
#endif // USE_MULTIPRECISION
