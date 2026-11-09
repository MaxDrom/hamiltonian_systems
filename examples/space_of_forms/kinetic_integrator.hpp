#pragma once
#include "auto_diff_base.hpp"
#include <auto_diff.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <precision.hpp>
using namespace boost::numeric::ublas;

inline auto get_kinetic() {
    using namespace AutoDiff;
    auto q1 = Variable<0>();
    auto q2 = Variable<1>();
    auto q3 = Variable<2>();
    auto p1 = Variable<3>();
    auto p2 = Variable<4>();
    auto p3 = Variable<5>();
    return Integer<2>() * sqrt(q1 * q1 + q2 * q2 + q3 * q3) *
           (p1 * p1 + p2 * p2 + p3 * p3);
}

inline void kinetic_solver(Real t, vector<Real> &x, Real h) {
    Real x0[3] = {x[0], x[1], x[2]};
    Real p0[3] = {x[3], x[4], x[5]};

    Real a = sqrt(x0[0] * x0[0] + x0[1] * x0[1] + x0[2] * x0[2]);
    if (a == 0)
        return;

    Real b = sqrt(p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2]);
    Real c = x0[0] * p0[0] + x0[1] * p0[1] + x0[2] * p0[2];

    Real Lx = x0[1] * p0[2] - x0[2] * p0[1];
    Real Ly = x0[2] * p0[0] - x0[0] * p0[2];
    Real Lz = x0[0] * p0[1] - x0[1] * p0[0];
    Real d = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

    Real tau = h;
    Real ab2 = a * b * b;
    Real H = 2 * ab2;

    if (d == 0) {
        // x0 || p0.
        Real r = a + 4 * c * tau + 4 * ab2 * tau * tau;
        Real sx = r / a;
        Real sp = (c + H * tau) / (a * r);
        for (int i = 0; i < 3; ++i) {
            x[i] = sx * x0[i];
            x[3 + i] = sp * x0[i];
        }
        return;
    }

    Real u = 2 * ab2 * tau + c;
    Real r = (u * u + d * d) / ab2;
    Real phi = 2 * atan2(2 * ab2 * tau * d, d * d + u * c);
    Real s = sin(phi);
    Real co = cos(phi);

    Real A = (r / a) * co - (r * c / (a * d)) * s;
    Real B = (a * r / d) * s;
    Real C = (H * tau / (a * r)) * co -
             ((a * a * b * b + c * H * tau) / (a * d * r)) * s;
    Real D = (a / r) * co + (a * (c + H * tau) / (r * d)) * s;

    for (int i = 0; i < 3; ++i) {
        x[i] = A * x0[i] + B * p0[i];
        x[3 + i] = C * x0[i] + D * p0[i];
    }
}
