#pragma once
#include "ham.hpp"
#include <auto_diff.hpp>
#include <auto_diff_base.hpp>
#include <auto_diff_expressions.hpp>
#include <boost/mpl/size_fwd.hpp>
#include <utility>

struct Nil {};

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
    using type = Concat_t<MakePairsWith_t<Head, Tail>,
                          typename MakeAllPairs<Tail>::type>;
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

template <typename List, size_t N_Bodies, size_t Dim> struct MakeKineticPart {};

template <size_t N_Bodies, size_t Dim>
struct MakeKineticPart<Nil, N_Bodies, Dim> {
    using type = AutoDiff::Constant<0.0>;
};

template <size_t Idx, double Mass, typename Tail, size_t N_Bodies, size_t Dim>
struct MakeKineticPart<List<BodyEntry<Idx, Mass>, Tail>, N_Bodies, Dim> {
    template <size_t... Ds>
    static decltype(auto) _sum_sq(std::index_sequence<Ds...>) {
        using namespace AutoDiff;
        return ((Variable<N_Bodies * Dim + Dim * Idx + Ds>{} *
                 Variable<N_Bodies * Dim + Dim * Idx + Ds>{}) +
                ...);
    }

    static decltype(auto) _make_part() {
        using namespace AutoDiff;
        return (_sum_sq(std::make_index_sequence<Dim>())) *
               Constant<1.0 / (2 * Mass)>{};
    }
    using _part_type = decltype(_make_part());
    using type =
        AutoDiff::Add<_part_type,
                      typename MakeKineticPart<Tail, N_Bodies, Dim>::type>;
};

template <typename List, size_t N_Bodies, size_t Dim>
struct MakePotentialPart {};

template <size_t N_Bodies, size_t Dim>
struct MakePotentialPart<Nil, N_Bodies, Dim> {
    using type = AutoDiff::Constant<0.0>;
};

template <size_t i, size_t j, double Mass_i, double Mass_j, typename Tail,
          size_t N_Bodies, size_t Dim>
struct MakePotentialPart<
    List<PairEntry<BodyEntry<i, Mass_i>, BodyEntry<j, Mass_j>>, Tail>, N_Bodies,
    Dim> {
    template <size_t... Ds>
    static decltype(auto) _sum_sq(std::index_sequence<Ds...>) {
        using namespace AutoDiff;
        return (((Variable<Dim * i + Ds>{} - Variable<Dim * j + Ds>{}) *
                 (Variable<Dim * i + Ds>{} - Variable<Dim * j + Ds>{})) +
                ...);
    }
    static decltype(auto) _make_part() {
        using namespace AutoDiff;
        auto r = sqrt(_sum_sq(std::make_index_sequence<Dim>()));
        return Constant<-1.0>{} / r;
    }
    using _part_type = decltype(_make_part());
    using type =
        AutoDiff::Add<_part_type,
                      typename MakePotentialPart<Tail, N_Bodies, Dim>::type>;
};

template <size_t Dim, double... Masses> struct NBodyProblem {
    static constexpr size_t N_Bodies = sizeof...(Masses);
    static constexpr size_t NCoords = 2 * N_Bodies * Dim;
    static constexpr size_t NFreedom = N_Bodies * Dim;
    decltype(auto) hamiltonian() const {
        using Bodies = MakeBodyList_t<Body<Masses>...>;
        using KineticPart =
            typename MakeKineticPart<Bodies, N_Bodies, Dim>::type;
        using Pairs = MakeAllPairs_t<Bodies>;
        using PotentialPart =
            typename MakePotentialPart<Pairs, N_Bodies, Dim>::type;
        return KineticPart{} + PotentialPart{};
    }

    decltype(auto) compiled_system() const {
        return Hamiltonian::compile_system<NFreedom>(hamiltonian());
    }
};

template <double... Masses>
decltype(auto) make_planar_nbody(Body<Masses>... args) {
    return NBodyProblem<2, Masses...>{};
}

template <double... Masses>
decltype(auto) make_general_nbody(Body<Masses>... args) {
    return NBodyProblem<3, Masses...>{};
}
