#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <nbody.hpp>

Real C(Real z) {
    if (z > 0)
        return (1 - cos(sqrt(z))) / z;
    if (z == 0)
        return Real(0.5);

    return (1 - cosh(sqrt(-z))) / z;
}

Real S(Real z) {
    if (z == 0)
        return Real(1) / Real(6);
    Real sq;
    if (z > 0) {
        sq = sqrt(z);

        return (sq - sin(sq)) / (sq * z);
    }

    sq = sqrt(-z);
    return (1 - sinh(sq)) / (sq * z);
}
using namespace boost::numeric::ublas;
using vec = vector_range<vector<Real>>;
void solveKeplerProblem(vec &r0, vec &v0, vec &r1, vec &v1, Real t) {
    auto v = v1 - v0;
    auto r = r1 - r0;
    Real rr = norm_2(r);
    Real vv = norm_2(v);
    // TODO: дописать аналитический решатель задачи кеплера
    r0 += v0 * t;
}
