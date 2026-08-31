#pragma once
#include <auto_diff_base.hpp>
#include <auto_diff_expressions.hpp>
#include <auto_diff_normalizer.hpp>
#include <auto_diff_squise.hpp>
#include <auto_diff_subst.hpp>

namespace AutoDiff {

template <IsExpression A, IsVariable Var> struct Diff {};

template <IsExpression A, IsExpression B, IsVariable Var>
struct Diff<Mul<A, B>, Var> {
    using _B = typename Diff<B, Var>::type;
    using _A = typename Diff<A, Var>::type;

    using type = Add<Mul<A, _B>, Mul<_A, B>>;
};

template <IsExpression A, IsExpression B, IsVariable Var>
struct Diff<Div<A, B>, Var> {
    using _B = typename Diff<B, Var>::type;
    using _A = typename Diff<A, Var>::type;

    using type =
        Mul<Minus<Mul<_A, B>, Mul<A, _B>>, Div<Constant<1.0>, Mul<B, B>>>;
};

template <double CC, size_t Idx> struct Diff<Constant<CC>, Variable<Idx>> {
    using type = Integer<0>;
};

template <int CC, size_t Idx> struct Diff<Integer<CC>, Variable<Idx>> {
    using type = Integer<0>;
};

template <IsExpression A, IsExpression B, IsVariable Var>
struct Diff<Add<A, B>, Var> {
    using _B = typename Diff<B, Var>::type;
    using _A = typename Diff<A, Var>::type;
    using type = Add<_A, _B>;
};

template <IsExpression A, IsExpression B, IsVariable Var>
struct Diff<Minus<A, B>, Var> {
    using _B = typename Diff<B, Var>::type;
    using _A = typename Diff<A, Var>::type;
    using type = Minus<_A, _B>;
};

template <size_t Idx> struct Diff<Variable<Idx>, Variable<Idx>> {
    using type = Integer<1>;
};

template <size_t Idx, size_t Idx2>
    requires(Idx != Idx2)
struct Diff<Variable<Idx>, Variable<Idx2>> {
  public:
    using type = Integer<0>;
};

template <IsExpression A, auto Func, size_t Idx>
struct DefaultFuncDiff : Expression<DefaultFuncDiff<A, Func, Idx>> {
    friend struct Expression<DefaultFuncDiff<A, Func, Idx>>;

  private:
    [[no_unique_address]]
    A _value;

    template <IsReal... Args> inline Real apply_impl(Args &&...args) const {
        Real h = std::sqrt(std::numeric_limits<Real>::epsilon());

        auto args_tuple = std::make_tuple(args...);

        auto args_plus = args_tuple;
        std::get<Idx>(args_plus) += h;

        auto args_minus = args_tuple;
        std::get<Idx>(args_minus) -= h;
        auto left = std::apply(_value, std::move(args_plus));
        auto right = std::apply(_value, std::move(args_minus));

        return (Func(left) - Func(right)) / (2 * h);
    }
    inline std::string print_impl() const {
        return "d_" + FuntionTagFromFunc<Func>{}.get_string();
    }
};

template <IsExpression A, auto Func> struct AnalyticalDiff {};

template <typename A, typename Var, auto Func>
concept HasAnalyticalDiff =
    requires { typename AnalyticalDiff<A, Func>::type; };

template <IsExpression A>
struct AnalyticalDiff<A, static_cast<Real (*)(Real &&)>(mp_sqrt)> {
    using type = Mul<Constant<0.5>, Div<Constant<1.0>, Sqrt<A>>>;
};

template <IsExpression A>
struct AnalyticalDiff<A, static_cast<Real (*)(Real &&)>(mp_cos)> {
    using type = Mul<Integer<-1>, Sin<A>>;
};

template <IsExpression A>
struct AnalyticalDiff<A, static_cast<Real (*)(Real &&)>(mp_sin)> {
    using type = Cos<A>;
};

template <IsExpression A>
struct AnalyticalDiff<A, static_cast<Real (*)(Real &&)>(mp_exp)> {
    using type = Exp<A>;
};

template <IsExpression A, auto Func, size_t Idx>
    requires(!HasAnalyticalDiff<A, Variable<Idx>, Func>)
struct Diff<Function<A, FuntionTagFromFunc<Func>>, Variable<Idx>> {
    using _A = typename Diff<A, Variable<Idx>>::type;
    using type = Mul<_A, DefaultFuncDiff<A, Func, Idx>>;
};

template <IsExpression A, auto Func, size_t Idx>
struct Diff<Function<A, FuntionTagFromFunc<Func>>, Variable<Idx>> {
    using _A = typename Diff<A, Variable<Idx>>::type;
    using _D = typename AnalyticalDiff<A, Func>::type;
    using type = Mul<_A, _D>;
};

template <IsExpression A, IsVariable Var>
FixedPointSquise<typename Diff<A, Var>::type>::type make_diff(const A &a,
                                                              const Var &b) {
    return {};
};
} // namespace AutoDiff
