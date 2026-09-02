#pragma once
#include "leapfrog.hpp"
#include "precision.hpp"
#include <utility>

namespace Yoshida {
using namespace boost::numeric::ublas;
namespace detail {

    template <size_t order> class YoshidaCoeffs {
      public:
        inline static const std::tuple<Real, Real> values = []() {
            Real x1 = Real(1.0) / (2 - mp_pow(2, Real(1) / Real(order - 1)));
            Real x0 = Real(1) - 2 * x1;

            for (size_t i = 0; i < 10; i++) {
                Real m11 = Real(1);
                Real m12 = Real(2);
                Real m21 = (order - 1) * mp_pow(x0, (order - 2));
                Real m22 = 2 * (order - 1) * mp_pow(x1, (order - 2));
                Real det = m11 * m22 - m12 * m21;
                Real f1 = x0 + 2 * x1 - 1;
                Real f2 = mp_pow(x0, (order - 1)) + 2 * mp_pow(x1, (order - 1));
                x0 -= (m22 * f1 - m12 * f2) / det;
                x1 -= (-m21 * f1 + m11 * f2) / det;
            }
            return std::make_tuple(x0, x1);
        }();
    };

    template <size_t order, typename h1_solver, typename h2_solver>
        requires IsSolverFunctor<h1_solver> && IsSolverFunctor<h2_solver>
    struct YoshidaImpl {
        using _LeapFrog = LeapFrog<h1_solver, h2_solver>;

        Real Step(Real t, vector<Real> &x, Real h, const _LeapFrog &leapfrog) {
            auto values = YoshidaCoeffs<order>::values;
            auto x0 = std::get<0>(values);
            auto x1 = std::get<1>(values);
            auto h0 = x0 * h;
            auto h1 = x1 * h;
            auto t_mid = _previos.Step(t, x, h1, leapfrog);
            t_mid = _previos.Step(t_mid, x, h0, leapfrog);
            return _previos.Step(t_mid, x, h1, leapfrog);
        }

      private:
        [[no_unique_address]]
        YoshidaImpl<order - 2, h1_solver, h2_solver> _previos;
    };

    template <typename h1_solver, typename h2_solver>
    struct YoshidaImpl<2, h1_solver, h2_solver> {

        using _LeapFrog = LeapFrog<h1_solver, h2_solver>;

        Real Step(Real t, vector<Real> &x, Real h, const _LeapFrog &leapfrog) {
            leapfrog.Step(t, x, h);
            return t + h;
        }
    };
} // namespace detail

template <size_t order, typename h1_solver, typename h2_solver>
    requires IsSolverFunctor<h1_solver> && IsSolverFunctor<h2_solver>
class Yoshida {

    using _LeapFrog = LeapFrog<h1_solver, h2_solver>;

  public:
    Yoshida(_LeapFrog &&leapfrog)
        : _leapfrog{std::forward<_LeapFrog>(leapfrog)} {}
    void Step(Real t, vector<Real> &x, Real h) {
        _impl.Step(t, x, h, _leapfrog);
    }

  private:
    _LeapFrog _leapfrog;
    detail::YoshidaImpl<order, h1_solver, h2_solver> _impl;
};
template <size_t order, typename h1_solver, typename h2_solver>
auto make_yoshida_from_leapfrog(LeapFrog<h1_solver, h2_solver> leapfrog) {
    return Yoshida<order, h1_solver, h2_solver>(
        std::forward<LeapFrog<h1_solver, h2_solver>>(leapfrog));
};
} // namespace Yoshida
