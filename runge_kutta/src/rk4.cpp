#include "runge_kutta.hpp"

namespace RungeKutta
{
    void RK4::Init(std::function<void(Real, const vector<Real>&, vector<Real>&)> rhs, int n)
    {
        _rhs = rhs;
        _n = n;
        k1 = vector<Real>(n);
        k2 = vector<Real>(n);
        k3 = vector<Real>(n);
        k4 = vector<Real>(n);
    }
    void RK4::Step(Real t, vector<Real>& x, Real h)
    {

        _rhs(t, x, k1);
        _rhs(t+h/2, x+h/2*k1, k2);
        _rhs(t+h/2, x+h/2*k2, k3);
        _rhs(t+h, x+h*k3, k4);

        x+= h*(k1+2*(k2+k3)+k4)/6;
    }
}