#pragma once
#include <boost/numeric/ublas/vector.hpp>
#include <precision.hpp>
namespace RungeKutta {
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;
template <typename Func>
concept IsRhs =
    requires(Func f, Real t, const vector<Real> &x, vector<Real> &result) {
        { f(t, x, result) } -> std::same_as<void>;
    };
// namespace detail

template <typename Solver>
concept IODESolver = requires(Solver solver, Real t, vector<Real> &x, Real h) {
    { solver.Step(t, x, h) } -> std::same_as<void>;
};

template <typename S>
concept HasReset = requires(S solver) {
    { solver.Reset() } -> std::same_as<void>;
};

template <typename Solver>
    requires IODESolver<Solver>
void Integrate(vector<Real> &x, Real t0, Real t1, Real h, Solver &integrator) {
    if constexpr (HasReset<Solver>)
        integrator.Reset();
    for (Real t = t0; t < t1; t += h) {
        integrator.Step(t, x, h);
    }
};

template <typename Solver>
    requires IODESolver<Solver>
void Integrate(vector<Real> &x, Real t0, Real t1, Real h, Solver &integrator,
               std::function<void(Real, const vector<Real> &)> callback) {
    if constexpr (HasReset<Solver>)
        integrator.Reset();
    Real t = t0;
    for (t = t0; t < t1; t += h) {
        callback(t, x);
        integrator.Step(t, x, h);
    }
    callback(t, x);
};
} // namespace RungeKutta
