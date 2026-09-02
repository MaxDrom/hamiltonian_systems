#pragma once
#include "auto_diff_base.hpp"
#include "leapfrog.hpp"
#include <auto_diff.hpp>
#include <cstddef>
#include <ham.hpp>
#include <utility>

namespace Yoshida {
using namespace AutoDiff;
template <size_t N, typename H, size_t... Is>
    requires IsExpression<H>
inline decltype(auto) _make_potential_integrator(H ham,
                                                 std::index_sequence<Is...>) {

    using P_Dots = std::tuple<decltype(make_diff(-ham, Variable<Is>{}))...>;
    return [=](Real t, vector<Real> &x, Real h) -> void {
        ((x[N + Is] =
              x[N + Is] +
              h * Hamiltonian::call_with_vector_impl(
                      t, x, typename std::tuple_element<Is, P_Dots>::type{},
                      std::make_index_sequence<2 * N>())),
         ...);
    };
}
template <size_t N, typename H, size_t... Is>
    requires IsExpression<H>
inline decltype(auto) _make_kinetic_integrator(H ham,
                                               std::index_sequence<Is...>) {

    using Q_Dots = std::tuple<decltype(make_diff(ham, Variable<Is + N>{}))...>;
    return [=](Real t, vector<Real> &x, Real h) -> void {
        ((x[Is] = x[Is] +
                  h * Hamiltonian::call_with_vector_impl(
                          t, x, typename std::tuple_element<Is, Q_Dots>::type{},
                          std::make_index_sequence<2 * N>())),
         ...);
    };
}
template <size_t NFreedom, typename H>
    requires IsExpression<H>
decltype(auto) make_leapfrog_for_sep(const H &ham) {
    auto kinetic_part = _make_kinetic_integrator<NFreedom>(
        ham, std::make_index_sequence<NFreedom>());
    auto potential_part = _make_potential_integrator<NFreedom>(
        ham, std::make_index_sequence<NFreedom>());
    return make_leapfrog(kinetic_part, potential_part);
}

} // namespace Yoshida
