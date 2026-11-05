#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <functional>
#include <newton.hpp>
#include <precision.hpp>
using namespace boost::numeric::ublas;
namespace Newton {
void jacobi(matrix<Real> &jac, const vector<Real> &x,
            std::function<void(const vector<Real> &, vector<Real> &)> f,
            const vector<Real> &f_x) {
    Real h = sqrt(std::numeric_limits<Real>::epsilon());

    int N = x.size();
    for (size_t i = 0; i < N; i++) {
        auto x_copy = x;
        x_copy(i) += h;
        vector<Real> new_f(x.size());
        f(x_copy, new_f);
        for (size_t j = 0; j < x.size(); j++) {
            jac(j, i) = (new_f(j) - f_x(j));
        }
    }
    jac = jac / h;
}

int solve_linear(const matrix<Real> &A, const vector<Real> &b,
                 vector<Real> &x) {
    matrix<Real> LU = A;
    permutation_matrix pm(A.size1());
    auto res = lu_factorize(LU, pm);
    if (res != 0)
        return res;
    x = b;
    lu_substitute(LU, pm, x);
    return 0;
}

int solve_linear(const matrix<Real> &A, const matrix<Real> &b,
                 matrix<Real> &x) {
    matrix<Real> LU = A;
    permutation_matrix pm(A.size1());
    auto res = lu_factorize(LU, pm);
    if (res != 0)
        return res;
    x = b;
    lu_substitute(LU, pm, x);
    return 0;
}

int inverse(const matrix<Real> &A, matrix<Real> &inverse) {
    matrix<Real> id = identity_matrix<Real>(A.size1());
    return solve_linear(A, id, inverse);
}

int solve(vector<Real> &x,
          std::function<void(const vector<Real> &, vector<Real> &)> f,
          Real damping) {
    vector<Real> dx(x.size());
    matrix<Real> jac(x.size(), x.size());
    vector<Real> fx(x.size());
    Real h = 10 * sqrt(std::numeric_limits<Real>::epsilon());
    do {
        f(x, fx);
        jacobi(jac, x, f, fx);
        int res = solve_linear(jac, fx, dx);
        if (res != 0) {
            if (norm_2(fx) <= std::numeric_limits<Real>::min())
                return 0;
            else
                return 1;
        }
        x -= damping * dx;

    } while (norm_2(dx) > 10 * std::max(h, norm_2(x) * h));

    return 0;
}

int solve(Real &x, std::function<Real(Real)> f) {
    Real h = 2 * sqrt(std::numeric_limits<Real>::epsilon());
    Real dx;
    Real fx;
    Real hh = sqrt(std::numeric_limits<Real>::epsilon());
    do {
        fx = f(x);
        auto df = (f(x + h) - fx) / h;
        if (mp_abs(df) < std::numeric_limits<Real>::min() || df == 0)
            return 1;
        dx = -fx / df;
        x = x + dx;

    } while (mp_abs(dx) > 10 * mp_max(hh, mp_abs(x) * hh));

    return 0;
}

int solve(Real &x, std::function<Real(Real)> f, std::function<Real(Real)> df) {

    Real dx;
    Real fx;
    do {
        fx = f(x);
        auto dfx = df(x);
        if (mp_abs(dfx) < std::numeric_limits<Real>::min() || dfx == 0)
            return 1;
        dx = -fx / dfx;
        x = x + dx;
        // std::cout<<x<<" "<<dx<<'\n';
    } while (mp_abs(fx) > std::numeric_limits<Real>::min() &&
             mp_abs(dx) > std::numeric_limits<Real>::epsilon());

    return 0;
}

Real find_minimum(Real a, Real b, std::function<Real(Real)> f) {
    Real phi = (sqrt(Real(5)) - Real(1)) / 2;
    Real c, d, x;
    Real eps = std::numeric_limits<Real>::epsilon();
    do {
        c = b - (b - a) * phi;
        d = a + (b - a) * phi;
        if (f(c) < f(d))
            b = d;
        else
            a = c;

        x = (a + b) / 2;

    } while ((b - a) > mp_max(eps, mp_abs(x) * eps));

    return x;
}
} // namespace Newton
