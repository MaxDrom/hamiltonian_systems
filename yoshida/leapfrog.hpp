#pragma once
#include <boost/numeric/ublas/vector.hpp>
#include <precision.hpp>
#include <runge_kutta/odesolver.hpp>
#include <type_traits>
#include <utility>
namespace Yoshida {
using namespace boost::numeric::ublas;
template <typename h1_solver, typename h2_sovler>
    requires RungeKutta::IsRhs<h1_solver> && RungeKutta::IsRhs<h2_sovler>
class LeapFrog {
    LeapFrog(h1_solver h1, h2_sovler h2)
        : _h1{std::forward<h1_solver>(h1)}, _h2{std::forward<h2_sovler>(h2)} {}
    void Step(Real t, vector<Real> &x) {
        _h1(t * Real(0.5), x);
        _h2(t, x);
        _h1(t * Real(0.5), x);
    }

  private:
    using CleahH1 = std::unwrap_ref_decay_t<h1_solver>;
    CleahH1 _h1;
    using CleahH2 = std::unwrap_ref_decay_t<h2_sovler>;
    CleahH2 _h2;
};
} // namespace Yoshida
