#include <auto_diff_base.hpp>
#include <auto_diff_expressions.hpp>
template <IsExpression Expr, IsExpression Pattern, IsExpression Target>
struct SubstWorker {
    using type = Expr;
};

template <IsExpression Pattern, IsExpression Target>
struct SubstWorker<Pattern, Pattern, Target> {
    using type = Target;
};

template <IsExpression A, IsExpression Pattern, IsExpression Target>
using Subst = ExpressionIterator<SubstWorker, A, Pattern, Target>;

template <IsExpression Expr, IsExpression Pattern, IsExpression Target>
typename Subst<Expr, Pattern, Target>::type
subst(const Expr &expr, const Pattern &pattern, const Target &target) {
    return {};
}