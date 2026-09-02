#pragma once
#include <boost/multiprecision/mpfr.hpp>
#include <cmath>
using Real = double;
using namespace boost::multiprecision;
inline Real mp_sin(Real value) { return std::sin(value); }

inline Real mp_cos(Real value) { return std::cos(value); }

inline Real mp_exp(Real value) { return std::exp(value); }

inline Real mp_sqrt(Real value) { return std::sqrt(value); }

inline Real mp_pow(Real value, Real deg) { return std::pow(value, deg); }
