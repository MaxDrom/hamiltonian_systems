#include "runge_kutta.hpp"

namespace RungeKutta
{
    void Integrate(int n, vector<Real>& x, Real t0, Real t1, Real h, 
                   std::function<void(Real, const vector<Real> &, vector<Real> &)> rhs, 
                   ODESolver &integrator)
    {
        integrator.Init(rhs, n);

        for(Real t = t0; t<t1; t+=h)
        {
            integrator.Step(t, x, h);
        }

    }
    void Integrate(int n, vector<Real>& x, Real t0, Real t1, Real h,
                   std::function<void(Real, const vector<Real> &, vector<Real> &)> rhs,
                   ODESolver &integrator,
                   std::function<void(Real, const vector<Real> &)> callback)
    {
        Real t = t0;
        integrator.Init(rhs, n);
        for(t = t0; t<t1; t+=h)
        {
            callback(t, x);
            integrator.Step(t, x, h);
        }
        callback(t, x);
    }
}