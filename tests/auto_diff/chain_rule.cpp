#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <auto_diff.hpp>
#include <auto_diff_base.hpp>
#include <auto_diff_expressions.hpp>
#include <doctest/doctest.h>
#include <precision.hpp>
using namespace AutoDiff;
TEST_CASE("Chain rule: d/dx sin(x^2)") {
    auto x = Variable<0>{};
    auto expr = sin(x * x);
    auto dexpr = make_diff(expr, x);
    for (int i = 0; i < 10; i++) {
        Real x0 = 2 * Real(3.14) / 10 * i;
        Real x2 = x0 * x0;
        Real expected = mp_cos(std::move(x2)) * 2 * x0;
        CHECK(dexpr(x0) == doctest::Approx(expected));
    }
}

TEST_CASE("Nested chain rule: d/dx exp(sin(x^2))") {
    auto x = Variable<0>{};
    auto expr = exp(sin(x * x));
    auto dexpr = make_diff(expr, x);
    Real x0 = 0.7;
    Real inner = x0 * x0;
    Real expected = mp_exp(std::move(mp_sin(std::move(inner)))) *
                    mp_cos(std::move(inner)) * 2 * x0;
    CHECK(dexpr(x0) == doctest::Approx(expected));
}

TEST_CASE("Chain rule through division and multiplication combined") {
    auto x = Variable<0>{};
    auto expr = Constant<1.0>{} / (Constant<1.0>{} + x * x);
    auto dexpr = make_diff(expr, x);
    Real x0 = 2.0;
    Real denom = 1 + x0 * x0;
    Real expected = -2 * x0 / (denom * denom);
    CHECK(dexpr(x0) == doctest::Approx(expected));
}
