#pragma once
#include "newton.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
#include <precision.hpp>
#include <string>
#include <type_traits>
namespace RungeKutta {
using namespace boost::multiprecision;
using namespace boost::numeric::ublas;
namespace detail {
    inline Real Legandre(Real x, int k) {
        if (k == 0)
            return 1;

        Real res1 = 1;
        auto res = x;
        for (int i = 1; i < k; i++) {
            auto tmp = res;
            res = (2 * i + 1) / Real(i + 1) * x * res - i / Real(i + 1) * res1;
            res1 = tmp;
        }
        return res;
    }

    inline Real DLegandre(Real x, int k) {
        if (k == 0)
            return 1;

        return k / (1 - x * x) * (Legandre(x, k - 1) - x * Legandre(x, k));
    }

} // namespace detail
class ODESolver {
  public:
    virtual void Init() = 0;
    virtual void Step(Real t, vector<Real> &x, Real h) = 0;
};
template <size_t n, auto _rhs> class RK4 final : public ODESolver {
  public:
    RK4() {};
    void Init() override {
        _n = n;
        k1 = vector<Real>(n);
        k2 = vector<Real>(n);
        k3 = vector<Real>(n);
        k4 = vector<Real>(n);
    }
    void Step(Real t, vector<Real> &x, Real h) override {
        _rhs(t, x, k1);
        _rhs(t + h / 2, x + h / 2 * k1, k2);
        _rhs(t + h / 2, x + h / 2 * k2, k3);
        _rhs(t + h, x + h * k3, k4);
        x += h * (k1 + 2 * (k2 + k3) + k4) / 6;
    }

  private:
    int _n;
    vector<Real> k1;
    vector<Real> k2;
    vector<Real> k3;
    vector<Real> k4;
};
template <size_t n, auto _rhs> class BRK final : public ODESolver {
  public:
    BRK(std::string path_to_data, int s) {
        _s = s;
        std::ifstream f(path_to_data);
        f >> c;
        f >> A;
        f >> b;
        f >> L;
        _rk4 = RK4<n, _rhs>();
    }
    BRK(int k, int s) {
        _s = s;
        auto N = k;
        vector<Real> roots1(N);
        Real pi = boost::math::constants::pi<Real>();
        for (size_t i = 1; i <= roots1.size(); i++) {
            roots1(i - 1) = cos((4 * i - 1) * (pi / (4 * N + 2)));
            Newton::solve(
                roots1(i - 1), [&N](Real x) { return detail::Legandre(x, N); },
                [&N](Real x) { return detail::DLegandre(x, N); });

            roots1(i - 1) = (roots1(i - 1) + 1) / 2;
        }

        vector<Real> roots(N);
        for (size_t i = 0; i < roots.size(); i++) {
            roots(i) = roots1(roots.size() - i - 1);
        }
        matrix<Real> Cm(N, N);
        matrix<Real> Cm_inv(N, N);
        matrix<Real> C(N, N);
        vector<Real> q(N);

        for (size_t i = 0; i < N; i++) {
            q(i) = 1 / Real(i + 1);
            for (size_t j = 0; j < N; j++) {
                Cm(i, j) = pow(roots(i), j);
                C(i, j) = Cm(i, j) * roots(i) / Real(j + 1);
            }
        }

        Newton::inverse(Cm, Cm_inv);
        auto A = prod(C, Cm_inv);
        auto b = prod(trans(Cm_inv), q);

        matrix<Real> L(N, N);
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                L(i, j) = detail::Legandre(2 * (1 + roots(i)) - 1, N) /
                          (2 * detail::DLegandre(2 * roots(j) - 1, N)) /
                          (1 + roots(i) - roots(j));
            }
        }
        c = roots;
        this->A = trans(A);
        this->b = b;
        this->L = L;
    }
    void Init() override {
        _first_step = true;
        K = matrix<Real>(n, c.size());
        y = matrix<Real>(n, c.size());
        buf = vector<Real>(n);
        buf2 = vector<Real>(n);
    }
    void Step(Real t, vector<Real> &x, Real h) override {
        {
            if (_first_step) {
                Prepare(t, x, h);
                _first_step = false;
            } else {
                y = h * prod(K, L);
                for (size_t i = 0; i < y.size2(); i++) {
                    column(y, i) += x;
                }
            }

            for (size_t i = 0; i < y.size2(); i++) {
                for (size_t j = 0; j < n; j++) {
                    buf(j) = y(j, i);
                    buf2(j) = K(j, i);
                }
                _rhs(t + c(i) * h, buf, buf2);
                for (size_t j = 0; j < n; j++) {
                    K(j, i) = buf2(j);
                }
            }

            for (int ss = 0; ss < _s; ss++) {

                y = h * prod(K, A);
                for (size_t i = 0; i < y.size2(); i++) {
                    column(y, i) += x;
                }

                for (size_t i = 0; i < y.size2(); i++) {
                    for (size_t j = 0; j < n; j++) {
                        buf(j) = y(j, i);
                        buf2(j) = K(j, i);
                    }
                    _rhs(t + c(i) * h, buf, buf2);
                    for (size_t j = 0; j < n; j++) {
                        K(j, i) = buf2(j);
                    }
                }
            }

            x = x + h * prod(K, b);
        }
    }

  private:
    void Prepare(Real t, const vector<Real> &x, Real h) {
        y = matrix<Real>(n, c.size());

        for (size_t i = 0; i < c.size(); i++) {
            _rk4.Init();

            buf = x;

            _rk4.Step(t, buf, c(i) * h);
            column(y, i) = buf;
        }
    }
    matrix<Real> y;
    matrix<Real> K;
    matrix<Real> A;
    matrix<Real> L;
    vector<Real> b;
    vector<Real> c;
    RK4<n, _rhs> _rk4;
    int _s;
    bool _first_step;
    vector<Real> buf;
    vector<Real> buf2;
};

template <typename Solver>
    requires std::is_base_of_v<ODESolver, Solver>
void Integrate(vector<Real> &x, Real t0, Real t1, Real h, Solver &integrator) {
    integrator.Init();

    for (Real t = t0; t < t1; t += h) {
        integrator.Step(t, x, h);
    }
};

template <typename Solver>
    requires std::is_base_of_v<ODESolver, Solver>
void Integrate(vector<Real> &x, Real t0, Real t1, Real h, Solver &integrator,
               std::function<void(Real, const vector<Real> &)> callback) {
    Real t = t0;
    integrator.Init();
    for (t = t0; t < t1; t += h) {
        callback(t, x);
        integrator.Step(t, x, h);
    }
    callback(t, x);
};
} // namespace RungeKutta
