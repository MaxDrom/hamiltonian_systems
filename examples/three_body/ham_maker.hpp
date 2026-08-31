#pragma once
#include <auto_diff.hpp>
#include <auto_diff_base.hpp>
#include <auto_diff_expressions.hpp>

template <typename Head, typename Tail> struct List {};

template <double mass> struct Body {};

template <size_t Idx, double Mass> struct BodyEntry {};

template <typename T> struct ExtractMass;

template <double Mass> struct ExtractMass<Body<Mass>> {
  static constexpr double value = Mass;
};

template <typename... Args> struct MakeList {
  using type = Nil;
};

template <typename Head, typename... Tail> struct MakeList<Head, Tail...> {
  using type = List<Head, typename MakeList<Tail...>::type>;
};

template <typename L1, typename L2> struct Concat;

template <typename L2> struct Concat<Nil, L2> {
  using type = L2;
};

template <typename Head, typename Tail, typename L2>
struct Concat<List<Head, Tail>, L2> {
  using type = List<Head, typename Concat<Tail, L2>::type>;
};

template <typename L1, typename L2>
using Concat_t = typename Concat<L1, L2>::type;

template <typename E1, typename E2> struct PairEntry {};

template <typename Head, typename L> struct MakePairsWith {
  using type = Nil;
};

template <typename Head, typename H2, typename Tail>
struct MakePairsWith<Head, List<H2, Tail>> {
  using type =
      List<PairEntry<Head, H2>, typename MakePairsWith<Head, Tail>::type>;
};

template <typename Head, typename L>
using MakePairsWith_t = typename MakePairsWith<Head, L>::type;

template <typename L> struct MakeAllPairs {
  using type = Nil;
};

template <typename Head, typename Tail> struct MakeAllPairs<List<Head, Tail>> {
  using type =
      Concat_t<MakePairsWith_t<Head, Tail>, typename MakeAllPairs<Tail>::type>;
};

template <typename L> using MakeAllPairs_t = typename MakeAllPairs<L>::type;

template <typename... Args> struct MakeBodyList {
private:
  template <std::size_t... Is>
  static auto make(std::index_sequence<Is...>) ->
      typename MakeList<BodyEntry<Is, ExtractMass<Args>::value>...>::type;

public:
  using type = decltype(make(std::index_sequence_for<Args...>{}));
};

template <typename... Args>
using MakeBodyList_t = typename MakeBodyList<Args...>::type;

template <typename List, size_t N_Bodies> struct MakeKineticPart {};

template <size_t N_Bodies> struct MakeKineticPart<Nil, N_Bodies> {
  using type = Constant<0.0>;
};

template <size_t Idx, double Mass, typename Tail, size_t N_Bodies>
struct MakeKineticPart<List<BodyEntry<Idx, Mass>, Tail>, N_Bodies> {
  static decltype(auto) _make_part() {
    auto _px = Variable<2 * N_Bodies + 2 * Idx>{};
    auto _py = Variable<2 * N_Bodies + 2 * Idx + 1>{};
    return (_px * _px + _py * _py) * Constant<1.0 / (2 * Mass)>{};
  }
  using _part_type = decltype(_make_part());
  using type = Add<_part_type, typename MakeKineticPart<Tail, N_Bodies>::type>;
};

template <typename List, size_t N_Bodies> struct MakePotentialPart {};

template <size_t N_Bodies> struct MakePotentialPart<Nil, N_Bodies> {
  using type = Constant<0.0>;
};

template <size_t i, size_t j, double Mass_i, double Mass_j, typename Tail,
          size_t N_Bodies>
struct MakePotentialPart<
    List<PairEntry<BodyEntry<i, Mass_i>, BodyEntry<j, Mass_j>>, Tail>,
    N_Bodies> {
  static decltype(auto) _make_part() {
    auto _qix = Variable<2 * i>{};
    auto _qiy = Variable<2 * i + 1>{};
    auto _qjx = Variable<2 * j>{};
    auto _qjy = Variable<2 * j + 1>{};
    auto r =
        sqrt((_qiy - _qjy) * (_qiy - _qjy) + (_qix - _qjx) * (_qix - _qjx));
    return Constant<-1.0>{} / r;
  }
  using _part_type = decltype(_make_part());
  using type =
      Add<_part_type, typename MakePotentialPart<Tail, N_Bodies>::type>;
};

template <double... Masses> decltype(auto) make_ham(Body<Masses>... args) {
  constexpr size_t N_Bodies = sizeof...(args);
  using Bodies = MakeBodyList_t<Body<Masses>...>;
  using Pairs = MakeAllPairs_t<Bodies>;
  using KineticPart = typename MakeKineticPart<Bodies, N_Bodies>::type;
  using PotentialPart = typename MakePotentialPart<Pairs, N_Bodies>::type;

  return KineticPart{} + PotentialPart{};
}
