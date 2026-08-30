#pragma once
#include <auto_diff_base.hpp>
#include <auto_diff_expressions.hpp>

template<IsExpression A>
struct Squiser
{
    using type = A;
};

template<IsNotNumber A>
struct Squiser<Mul<Integer<1>, A>> { using type = A; };

template<IsNotNumber A>
struct Squiser<Mul<Constant<1.0>, A>> { using type = A; };

template<IsNotNumber A>
struct Squiser<Mul<A, Integer<1>>> { using type = A; };

template<IsNotNumber A>
struct Squiser<Mul<A, Constant<1.0>>> { using type = A; };


template<IsNotNumber A>
struct Squiser<Add<Integer<0>, A>> { using type = A; };

template<IsNotNumber A>
struct Squiser<Add<Constant<0.0>, A>> { using type = A; };

template<IsNotNumber A>
struct Squiser<Add<A, Integer<0>>> { using type = A; };

template<IsNotNumber A>
struct Squiser<Add<A, Constant<0.0>>> { using type = A; };

template<int L, int R>
struct Squiser<Mul<Integer<L>, Integer<R>>> { using type = Integer<L*R>; };

template<double L, double R>
struct Squiser<Mul<Constant<L>, Constant<R>>> { using type = Constant<L*R>; };

template<int L, double R>
struct Squiser<Mul<Integer<L>, Constant<R>>> { using type = Constant<L*R>; };
template<double L, int R>
struct Squiser<Mul<Constant<L>, Integer<R>>> { using type = Constant<L*R>; };

template<double L, double R>
struct Squiser<Div<Constant<L>, Constant<R>>> { using type = Constant<L/R>; };

template<int L, double R>
struct Squiser<Div<Integer<L>, Constant<R>>> { using type = Constant<L/R>; };
template<double L, int R>
struct Squiser<Div<Constant<L>, Integer<R>>> { using type = Constant<L/R>; };

template<int L, int R>
struct Squiser<Add<Integer<L>, Integer<R>>> { using type = Integer<L+R>; };
template<double L, double R>
struct Squiser<Add<Constant<L>, Constant<R>>> { using type = Constant<L+R>; };
template<int L, double R>
struct Squiser<Add<Integer<L>, Constant<R>>> { using type = Constant<L+R>; };
template<double L, int R>
struct Squiser<Add<Constant<L>, Integer<R>>> { using type = Constant<L+R>; };

template<int L, int R>
struct Squiser<Minus<Integer<L>, Integer<R>>> { using type = Integer<L-R>; };
template<double L, double R>
struct Squiser<Minus<Constant<L>, Constant<R>>> { using type = Constant<L-R>; };
template<int L, double R>
struct Squiser<Minus<Integer<L>, Constant<R>>> { using type = Constant<L-R>; };
template<double L, int R>
struct Squiser<Minus<Constant<L>, Integer<R>>> { using type = Constant<L-R>; };

template<IsNotNumber A>
struct Squiser<Mul<A, Div<Constant<1.0>, A>>>
{
    using type = Constant<1.0>;
};

template<IsNotNumber A>
struct Squiser<Mul<Div<Constant<1.0>, A>, A>>
{
    using type = Constant<1.0>;
};

template<IsNotNumber A>
struct Squiser<Mul<A, Div<Integer<1>, A>>>
{
    using type = Constant<1.0>;
};

template<IsNotNumber A>
struct Squiser<Mul<Div<Integer<1>, A>, A>>
{
    using type = Constant<1.0>;
};

template<IsNotNumber expr, int Val>
struct Squiser<Minus<expr, Integer<Val>>>
{
    using type = Add<Integer<-Val>, expr>;
};

template<IsNotNumber expr, double Val>
struct Squiser<Minus<expr, Constant<Val>>>
{
    using type = Add<Constant<-Val>, expr>;
};

template<IsNotNumber A>
struct Squiser<Add<A, A>> {
    using type = Mul<Integer<2>, A>;
};

template<IsNotNumber A>
struct Squiser<Mul<Constant<0.0>, A>> {
    using type = Integer<0>;
};

template<IsNotNumber A>
struct Squiser<Mul<Integer<0>, A>> {
    using type = Integer<0>;
};

template<IsNotNumber A>
struct Squiser<Mul<A, Constant<0.0>>> {
    using type = Integer<0>;
};

template<IsNotNumber A>
struct Squiser<Mul<A, Integer<0>>> {
    using type = Integer<0>;
};

template<IsNotNumber A, size_t Val>
struct Squiser<Add<Mul<Integer<Val>, A>, A>>
{
    using type = Mul<Integer<Val+1>, A>;
};

template<IsNotNumber A, size_t Val>
struct Squiser<Add<A, Mul<Integer<Val>, A>>>
{
    using type = Mul<Integer<Val+1>, A>;
};

template<IsNotNumber A, size_t Val>
struct Squiser<Add<A, Mul<A, Integer<Val>>>>
{
    using type = Mul<Integer<Val+1>, A>;
};

template<IsNotNumber A, size_t Val>
struct Squiser<Add<Mul<A, Integer<Val>>, A>>
{
    using type = Mul<Integer<Val+1>, A>;
};

template<IsNotNumber A>
struct Squiser<Minus<Integer<0>, A>>
{
    using type = Mul<Integer<-1>, A>;
};

template<IsNotNumber A>
struct Squiser<Minus<Constant<0.0>, A>>
{
    using type = Mul<Integer<-1>, A>;
};

template<IsNotNumber A>
struct Squiser<Minus<A, Integer<0>>>
{
    using type = A;
};

template<IsNotNumber A>
struct Squiser<Minus<A, Constant<0.0>>>
{
    using type = A;
};

template<IsExpression A>
struct Squiser<Mul<Sqrt<A>, Sqrt<A>>>
{
    using type = A;
};

template<IsExpression A>
struct Squiser<Mul<Sqrt<A>, Div<Constant<1.0>, A>>>
{
    using type = Div<Constant<1.0>,Sqrt<A>>;
};

