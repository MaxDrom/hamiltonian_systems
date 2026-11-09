#include "auto_diff_base.hpp"
#include <auto_diff.hpp>
using namespace AutoDiff;

template <double m1, double m2, double m3> decltype(auto) get_potential() {
    constexpr double M = m1 + m2 + m3;
    auto M2 = sqrt(Constant<m1 * m2 * m3 * M>());
    auto q0 = Variable<0>{};
    auto q1 = Variable<1>{};
    auto q2 = Variable<2>{};

    auto rho = sqrt(q0 * q0 + q1 * q1 + q2 * q2);
    auto r12 = sqrt(Constant<(m1 + m2) / (2 * m1 * m2)>() * (rho + q0));
    auto r13 =
        sqrt(Constant<(m1 + m3) / (2 * m1 * m3)>() * rho +
             q0 * Constant<(m2 * m3 - m1 * M) / (2 * m1 * m3 * (M - m3))>() +
             q1 * M2 / Constant<(m1 * m3 * (M - m3))>());
    auto r23 =
        sqrt(Constant<(m2 + m3) / (2 * m2 * m3)>() * rho +
             q0 * Constant<(m1 * m3 - m2 * M) / (2 * m2 * m3 * (M - m3))>() -
             q1 * M2 / Constant<(m2 * m3 * (M - m3))>());

    auto V = Constant<m1 * m2>() / r12 + Constant<m2 * m3>() / r23 +
             Constant<m1 * m3>() / r13;

    return -V;
}
