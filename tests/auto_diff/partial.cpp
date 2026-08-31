#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <auto_diff.hpp>
#include <auto_diff_base.hpp>
#include <auto_diff_expressions.hpp>
#include <doctest/doctest.h>
#include <precision.hpp>

TEST_CASE("Partial derivative ignores other variables") {
  // f(x, y) = x^2 * y
  auto x = Variable<0>{};
  auto y = Variable<1>{};
  auto expr = x * x * y;

  auto df_dx = make_diff(expr, x);
  auto df_dy = make_diff(expr, y);

  Real x0 = 2.0, y0 = 3.0;
  CHECK(df_dx(x0, y0) == doctest::Approx(2 * x0 * y0));
  CHECK(df_dy(x0, y0) == doctest::Approx(x0 * x0));
}

TEST_CASE(
    "Mixed second partials are equal (Schwarz theorem) for smooth expression") {
  // f(x,y) = sin(x*y)
  auto x = Variable<0>{};
  auto y = Variable<1>{};
  auto expr = sin(x * y);

  auto d2f_dxdy = make_diff(make_diff(expr, x), y);
  auto d2f_dydx = make_diff(make_diff(expr, y), x);

  Real x0 = 1.1, y0 = 0.6;
  CHECK(d2f_dxdy(x0, y0) == doctest::Approx(d2f_dydx(x0, y0)));
}
