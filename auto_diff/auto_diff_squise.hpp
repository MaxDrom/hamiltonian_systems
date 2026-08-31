#pragma once
#include <auto_diff_base.hpp>
#include <auto_diff_expressions.hpp>
#include <auto_diff_normalizer.hpp>
#include <auto_diff_squiser_rules.hpp>
#include <type_traits>
namespace AutoDiff {

template <IsExpression A> using Squise = ExpressionIterator<Squiser, A>;

template <typename Current, typename Next, bool IsSame> struct FixedPointImpl {
    using type = Current;
};

template <typename Current, typename Next>
struct FixedPointImpl<Current, Next, false> {
    using _squized = typename Squise<typename Normalize<Next>::type>::type;
    using type = typename FixedPointImpl<Next, _squized,
                                         std::is_same_v<Next, _squized>>::type;
};

template <IsExpression T> struct FixedPointSquise {
    using FirstPass = typename Squise<typename Normalize<T>::type>::type;

    using type = typename FixedPointImpl<T, FirstPass,
                                         std::is_same_v<T, FirstPass>>::type;
};

template <IsExpression Expr>
typename FixedPointSquise<Expr>::type simplify(const Expr &expr) {
    return {};
}
} // namespace AutoDiff
