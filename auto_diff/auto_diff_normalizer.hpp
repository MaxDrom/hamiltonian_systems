#pragma once
#include "auto_diff_base.hpp"
#include "auto_diff_expressions.hpp"

template <IsBinaryOpTag OpTag> constexpr double combine(double a, double b) {
  return 0;
}

template <IsBinaryOpTag OpTag> constexpr double identity() { return 0.0; }

template <> constexpr double combine<AddTag>(double a, double b) {
  return a + b;
}

template <> constexpr double combine<MulTag>(double a, double b) {
  return a * b;
}

template <> constexpr double identity<AddTag>() { return 0.0; }

template <> constexpr double identity<MulTag>() { return 1.0; }

template <typename Head, typename Tail> struct TypeList {};

struct Nil {};

template <IsExpression A> struct Normalize {
  using type = A;
};

template <IsExpression A, IsFunctionTag FuncTag>
struct Normalize<Function<A, FuncTag>> {
  using _A = typename Normalize<A>::type;
  using type = Function<_A, FuncTag>;
};

template <IsExpression Left, IsExpression Right, IsBinaryOpTag OpTag>
struct Normalize<BinaryOp<Left, Right, OpTag>> {
  using _Left = typename Normalize<Left>::type;
  using _Right = typename Normalize<Right>::type;
  using type = BinaryOp<_Left, _Right, OpTag>;
};

struct EmptyAcc {};

template <IsBinaryOpTag op, typename Acc, typename List, double AccConstant>
struct FoldLeftList;

template <IsBinaryOpTag op, IsExpression Acc, double AccConstant>
struct FoldLeftList<op, Acc, Nil, AccConstant> {
  using type = BinaryOp<Acc, Constant<AccConstant>, op>;
};

template <IsBinaryOpTag op, double AccConstant>
struct FoldLeftList<op, EmptyAcc, Nil, AccConstant> {
  using type = Constant<AccConstant>;
};

template <IsBinaryOpTag op, IsNotNumber Head, typename Tail, double AccConstant>
struct FoldLeftList<op, EmptyAcc, TypeList<Head, Tail>, AccConstant> {
  using NewAcc = Head;
  using type = typename FoldLeftList<op, NewAcc, Tail, AccConstant>::type;
};

template <IsBinaryOpTag op, IsExpression Acc, IsNotNumber Head, typename Tail,
          double AccConstant>
struct FoldLeftList<op, Acc, TypeList<Head, Tail>, AccConstant> {
  using NewAcc = BinaryOp<Acc, Head, op>;
  using type = typename FoldLeftList<op, NewAcc, Tail, AccConstant>::type;
};

template <IsBinaryOpTag op, typename Acc, int Head, typename Tail,
          double AccConstant>
struct FoldLeftList<op, Acc, TypeList<Integer<Head>, Tail>, AccConstant> {
  using type =
      typename FoldLeftList<op, Acc, Tail,
                            combine<op>((double)Head, AccConstant)>::type;
};

template <IsBinaryOpTag op, typename Acc, double Head, typename Tail,
          double AccConstant>
struct FoldLeftList<op, Acc, TypeList<Constant<Head>, Tail>, AccConstant> {
  using type = typename FoldLeftList<op, Acc, Tail,
                                     combine<op>(Head, AccConstant)>::type;
};

template <IsBinaryOpTag op, typename List> struct FoldLeft {
  using type = typename FoldLeftList<op, EmptyAcc, List, identity<op>()>::type;
};

template <IsBinaryOpTag OpTag, IsExpression Node, typename Tail>
struct ToFlatten {
  using type = TypeList<typename Normalize<Node>::type, Tail>;
};

template <IsBinaryOpTag OpTag, IsExpression L, IsExpression R, typename Tail>
struct ToFlatten<OpTag, BinaryOp<L, R, OpTag>, Tail> {
  using type =
      ToFlatten<OpTag, L, typename ToFlatten<OpTag, R, Tail>::type>::type;
};

template <IsExpression L, IsExpression R>
struct Normalize<BinaryOp<L, R, AddTag>> {
  using type = typename FoldLeft<
      AddTag,
      typename ToFlatten<AddTag, BinaryOp<L, R, AddTag>, Nil>::type>::type;
};

template <IsExpression L, IsExpression R>
struct Normalize<BinaryOp<L, R, MulTag>> {
  using type = typename FoldLeft<
      MulTag,
      typename ToFlatten<MulTag, BinaryOp<L, R, MulTag>, Nil>::type>::type;
};
