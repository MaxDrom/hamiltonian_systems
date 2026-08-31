#pragma once
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <functional>
#include <precision.hpp>
#include <string>
#include <type_traits>
namespace RungeKutta {
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;

class ODESolver {
  public:
    virtual void
    Init(std::function<void(Real, const vector<Real> &, vector<Real> &)> rhs,
         int n) = 0;
    virtual void Step(Real t, vector<Real> &x, Real h) = 0;
};

class RK4 final : public ODESolver {
  public:
    RK4() {};
    void
    Init(std::function<void(Real, const vector<Real> &, vector<Real> &)> rhs,
         int n) override;
    void Step(Real t, vector<Real> &x, Real h) override;

  private:
    std::function<void(Real, const vector<Real> &, vector<Real> &)> _rhs;
    int _n;
    vector<Real> k1;
    vector<Real> k2;
    vector<Real> k3;
    vector<Real> k4;
};

class BRK final : public ODESolver {
  public:
    BRK(std::string path_to_data, int s);
    BRK(int k, int s);
    void
    Init(std::function<void(Real, const vector<Real> &, vector<Real> &)> rhs,
         int n) override;
    void Step(Real t, vector<Real> &x, Real h) override;

  private:
    void Prepare(Real t, const vector<Real> &x, Real h);
    std::function<void(Real, const vector<Real> &, vector<Real> &)> _rhs;
    int _n;
    matrix<Real> y;
    matrix<Real> K;
    matrix<Real> A;
    matrix<Real> L;
    vector<Real> b;
    vector<Real> c;
    RK4 _rk4;
    int _s;
    bool _first_step;
    vector<Real> buf;
    vector<Real> buf2;
};

template <typename Solver>
    requires std::is_base_of_v<ODESolver, Solver>
void Integrate(
    int n, vector<Real> &x, Real t0, Real t1, Real h,
    std::function<void(Real, const vector<Real> &, vector<Real> &)> rhs,
    Solver &integrator) {
    integrator.Init(rhs, n);

    for (Real t = t0; t < t1; t += h) {
        integrator.Step(t, x, h);
    }
};

template <typename Solver>
    requires std::is_base_of_v<ODESolver, Solver>
void Integrate(
    int n, vector<Real> &x, Real t0, Real t1, Real h,
    std::function<void(Real, const vector<Real> &, vector<Real> &)> rhs,
    Solver &integrator,
    std::function<void(Real, const vector<Real> &)> callback) {
    Real t = t0;
    integrator.Init(rhs, n);
    for (t = t0; t < t1; t += h) {
        callback(t, x);
        integrator.Step(t, x, h);
    }
    callback(t, x);
};
} // namespace RungeKutta
