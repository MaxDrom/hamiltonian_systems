#pragma once
#include <auto_diff.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace Hamiltonian {
using namespace boost::numeric::ublas;
template <typename Func, std::size_t... I>
    requires std::invocable<Func, decltype((I, Real{}))..., Real>
decltype(auto)
    call_with_vector_impl(Real t,
                          const boost::numeric::ublas::vector<Real> &vec,
                          Func &&f, std::index_sequence<I...>) {
    return std::invoke(std::forward<Func>(f), vec[I]..., t);
}

template <size_t N, IsExpression H, size_t... Is>
inline decltype(auto) _compile_system(H ham, std::index_sequence<Is...>) {

    using Q_Dots = std::tuple<decltype(make_diff(ham, Variable<Is + N>{}))...>;
    using P_Dots = std::tuple<decltype(make_diff(-ham, Variable<Is>{}))...>;
    return [=](Real t, const vector<Real> &x, vector<Real> &result) -> void {
        ((result[Is] = call_with_vector_impl(
              t, x, typename std::tuple_element<Is, Q_Dots>::type{},
              std::make_index_sequence<2 * N>())),
         ...);
        ((result[Is + N] = call_with_vector_impl(
              t, x, typename std::tuple_element<Is, P_Dots>::type{},
              std::make_index_sequence<2 * N>())),
         ...);
    };
}

template <size_t N, IsExpression H> decltype(auto) compile_system(H ham) {
    return _compile_system<N>(ham, std::make_index_sequence<N>());
}
} // namespace Hamiltonian
