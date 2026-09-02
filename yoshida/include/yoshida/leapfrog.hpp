#pragma once
#include <boost/numeric/ublas/vector.hpp>
#include <precision.hpp>
#include <runge_kutta/odesolver.hpp>
#include <type_traits>
#include <utility>
namespace Yoshida {
using namespace boost::numeric::ublas;

template <typename F>
concept IsSolverFunctor = requires(F f, Real t, vector<Real> &x, Real h) {
    { f(t, x, h) } -> std::same_as<void>;
};
template <typename h1_solver, typename h2_solver>
    requires IsSolverFunctor<h1_solver> && IsSolverFunctor<h2_solver>
struct LeapFrog {
    LeapFrog(h1_solver &&h1, h2_solver &&h2)
        : _h1{std::forward<h1_solver>(h1)}, _h2{std::forward<h2_solver>(h2)} {}
    void Step(Real t, vector<Real> &x, Real h) const {
        auto h_half = h * Real(0.5);
        auto t_mid = t + h_half;
        _h1(t, x, h_half);
        _h2(t_mid, x, h);
        _h1(t_mid, x, h_half);
    }

  private:
    using CleahH1 = std::unwrap_ref_decay_t<h1_solver>;
    CleahH1 _h1;
    using CleahH2 = std::unwrap_ref_decay_t<h2_solver>;
    CleahH2 _h2;
};

template <typename h1_solver, typename h2_solver>
    requires IsSolverFunctor<h1_solver> && IsSolverFunctor<h2_solver>
auto make_leapfrog(h1_solver &&h1, h2_solver &&h2) {
    return LeapFrog<h1_solver, h2_solver>{std::forward<h1_solver>(h1),
                                          std::forward<h2_solver>(h2)};
}

} // namespace Yoshida
