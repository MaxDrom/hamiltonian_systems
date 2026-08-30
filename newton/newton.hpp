#pragma once
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <functional>
#include <precision.hpp>
namespace Newton
{
    using namespace boost::numeric::ublas;
    /*
    @param x inital value
    @param f equation to solve
    */
    int solve(vector<Real> &x, std::function<void(const vector<Real> &, vector<Real>&)> f, Real damping = 1);
    int solve_linear(const matrix<Real> &A, const vector<Real> &b, vector<Real> &x);
    int solve_linear(const matrix<Real> &A, const matrix<Real> &b, matrix<Real> &x);
    int solve(Real &x, std::function<Real(Real)> f);
    int solve(Real &x, std::function<Real(Real)> f, std::function<Real(Real)> df);
    int inverse(const matrix<Real> &A, matrix<Real> &inverse);
    Real find_minimum(Real a, Real b, std::function<Real(Real)> f);
}