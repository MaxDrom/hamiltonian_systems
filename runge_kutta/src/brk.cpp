#include "runge_kutta.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <fstream>
#include <newton.hpp>

Real Legandre(Real x, int k) {
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

Real DLegandre(Real x, int k) {
    if (k == 0)
        return 1;

    return k / (1 - x * x) * (Legandre(x, k - 1) - x * Legandre(x, k));
}

namespace RungeKutta {
BRK::BRK(std::string path_to_data, int s) {
    _s = s;
    std::ifstream f(path_to_data);
    f >> c;
    f >> A;
    f >> b;
    f >> L;
    _rk4 = RK4();
}

BRK::BRK(int k, int s) {
    _s = s;
    auto N = k;
    vector<Real> roots1(N);
    Real pi = boost::math::constants::pi<Real>();
    for (size_t i = 1; i <= roots1.size(); i++) {
        roots1(i - 1) = cos((4 * i - 1) * (pi / (4 * N + 2)));
        Newton::solve(
            roots1(i - 1), [&N](Real x) { return Legandre(x, N); },
            [&N](Real x) { return DLegandre(x, N); });

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
            L(i, j) = Legandre(2 * (1 + roots(i)) - 1, N) /
                      (2 * DLegandre(2 * roots(j) - 1, N)) /
                      (1 + roots(i) - roots(j));
        }
    }
    c = roots;
    this->A = trans(A);
    this->b = b;
    this->L = L;
}

void BRK::Init(
    std::function<void(Real, const vector<Real> &, vector<Real> &)> rhs,
    int n) {
    _rhs = rhs;
    _first_step = true;
    _n = n;

    K = matrix<Real>(n, c.size());
    y = matrix<Real>(n, c.size());
    buf = vector<Real>(n);
    buf2 = vector<Real>(n);
}

void BRK::Step(Real t, vector<Real> &x, Real h) {
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
        for (size_t j = 0; j < _n; j++) {
            buf(j) = y(j, i);
            buf2(j) = K(j, i);
        }
        _rhs(t + c(i) * h, buf, buf2);
        for (size_t j = 0; j < _n; j++) {
            K(j, i) = buf2(j);
        }
    }

    for (int ss = 0; ss < _s; ss++) {

        y = h * prod(K, A);
        for (size_t i = 0; i < y.size2(); i++) {
            column(y, i) += x;
        }

        for (size_t i = 0; i < y.size2(); i++) {
            for (size_t j = 0; j < _n; j++) {
                buf(j) = y(j, i);
                buf2(j) = K(j, i);
            }
            _rhs(t + c(i) * h, buf, buf2);
            for (size_t j = 0; j < _n; j++) {
                K(j, i) = buf2(j);
            }
        }
    }

    x = x + h * prod(K, b);
}

void BRK::Prepare(Real t, const vector<Real> &x, Real h) {
    y = matrix<Real>(_n, c.size());

    for (size_t i = 0; i < c.size(); i++) {
        _rk4.Init(_rhs, _n);

        buf = x;

        _rk4.Step(t, buf, c(i) * h);
        column(y, i) = buf;
    }
}
} // namespace RungeKutta
