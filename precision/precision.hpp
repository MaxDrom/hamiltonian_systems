#pragma once
#include <boost/math/constants/constants.hpp>
#ifndef USE_MULTIPRECISION
using Real = double;
#endif

#ifdef USE_MULTIPRECISION

#include <boost/multiprecision/fwd.hpp>
#include <boost/multiprecision/mpfr.hpp>
using namespace boost::multiprecision;
using Real = boost::multiprecision::static_mpfr_float_50;

#endif // USE_MULTIPRECISION
