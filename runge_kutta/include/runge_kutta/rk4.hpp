#pragma once
#include "odesolver.hpp"
namespace RungeKutta {
template <size_t n, typename Rhs>
    requires IsRhs<Rhs>
class RK4 {
  public:
    RK4(Rhs &&rhs)
        : _rhs{std::forward<Rhs>(rhs)}, k1{n}, k2{n}, k3{n}, k4{n} {};
    void Step(Real t, vector<Real> &x, Real h) {
        _rhs(t, x, k1);
        _rhs(t + h / 2, x + h / 2 * k1, k2);
        _rhs(t + h / 2, x + h / 2 * k2, k3);
        _rhs(t + h, x + h * k3, k4);
        x += h * (k1 + 2 * (k2 + k3) + k4) / 6;
    }

  private:
    using CleanRhs = std::unwrap_ref_decay_t<Rhs>;
    CleanRhs _rhs;
    vector<Real> k1;
    vector<Real> k2;
    vector<Real> k3;
    vector<Real> k4;
};
template <size_t n, typename Rhs>
    requires IsRhs<Rhs>
auto make_rk4(Rhs &&rhs) {
    return RK4<n, Rhs>(std::forward<Rhs>(rhs));
}
} // namespace RungeKutta
