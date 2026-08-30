#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "auto_diff_base.hpp"
#include <auto_diff.hpp>
#include <doctest/doctest.h>
#include <precision.hpp>

TEST_CASE("Derivative of constant is zero") {
    auto expr = Constant<5.0>{};
    auto x = Variable<0>{};
    auto dexpr = make_diff(expr, x);
    CHECK(dexpr(Real(3.0)) == doctest::Approx(0.0));
}

TEST_CASE("Derivative of variable w.r.t. itself is one") {
    auto x = Variable<0>{};
    auto expr = x;
    auto dexpr =make_diff(expr, x);
    CHECK(dexpr(Real(3.0)) == doctest::Approx(1.0));
}

TEST_CASE("Derivative of variable w.r.t. different variable is zero") {
    auto x = Variable<0>{};
    auto y = Variable<1>{};
    auto expr = x;
    auto dexpr = make_diff(expr, y);
    CHECK(dexpr(Real(3.0), Real(2.0)) == doctest::Approx(0.0));
}

TEST_CASE("Sum rule: d/dx(x + x^2) = 1 + 2x") {
    auto x = Variable<0>{};
    auto expr = x + x*x;
    auto dexpr = make_diff(expr, x);
    CHECK(dexpr(Real(3.0)) == doctest::Approx(7.0)); // 1 + 2*3
}

TEST_CASE("Product rule: d/dx(x * sin(x))") {
    auto x = Variable<0>{};
    auto expr = x * sin(x);
    auto dexpr = make_diff(expr, x);
    Real x0 = 1.5;
    Real expected = mp_sin(std::move(x0)) + x0 * mp_cos(std::move(x0));
    CHECK(dexpr(x0) == doctest::Approx(expected));
}

TEST_CASE("Quotient rule: d/dx(sin(x)/x)") {
    auto x = Variable<0>{};
    auto expr = sin(x) / x;
    auto dexpr = make_diff(expr, x);
    Real x0 = 2.0;
    Real expected = (std::cos(x0)*x0 - std::sin(x0)) / (x0*x0);
    CHECK(dexpr(x0) == doctest::Approx(expected));
}

