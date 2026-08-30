#include <type_traits>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <auto_diff_base.hpp>
#include <auto_diff_expressions.hpp>
#include <auto_diff.hpp>
#include <doctest/doctest.h>
#include <precision.hpp>

template<unsigned long long Seed>
struct CTRand {
    static constexpr unsigned long long value = (Seed * 6364136223846793005ULL + 1442695040888963407ULL);
    static constexpr unsigned long long next_seed = value;
};

template<unsigned long long Seed, int Depth, int Choice>
struct RandomChoices {};

template<unsigned long long Seed, int Depth>
struct RandomExpr {
    static constexpr unsigned long long r = CTRand<Seed>::value;
    static constexpr int choice = (r >> 32) % 4;

    using type = std::conditional_t<Depth == 0,
        Variable<0>,
        typename RandomChoices<Seed, Depth, choice>::type>;
};

template<unsigned long long Seed, int Depth>
struct RandomChoices<Seed, Depth, 0> {
    static constexpr unsigned long long next_seed = CTRand<Seed>::next_seed;
    using type = Mul<typename RandomExpr<next_seed, Depth-1>::type,
                typename RandomExpr<CTRand<next_seed>::next_seed, Depth-1>::type>;
};

template<unsigned long long Seed, int Depth>
struct RandomChoices<Seed, Depth, 1> {
    using type = Sin<typename RandomExpr<CTRand<Seed>::next_seed, Depth-1>::type>;
};

template<unsigned long long Seed, int Depth>
struct RandomChoices<Seed, Depth, 2> {
    using type = Cos<typename RandomExpr<CTRand<Seed>::next_seed, Depth-1>::type>;
};

template<unsigned long long Seed, int Depth>
struct RandomChoices<Seed, Depth, 3> {
    using type = Exp<typename RandomExpr<CTRand<Seed>::next_seed, Depth-1>::type>;
};

template<unsigned long long Seed>
struct RandomExpr<Seed, 0> {
    using type = Variable<0>;
};

template<IsExpression A>
inline Real central_difference(const A& expr, Real x, Real h = 1e-6)
{
    return (expr(x+h) - expr(x-h))/(2*h);
}

TEST_CASE_TEMPLATE("Random expression from seed", Seed,
                    std::integral_constant<unsigned long long, 1>,
                    std::integral_constant<unsigned long long, 42>,
                    std::integral_constant<unsigned long long, 145>,
                    std::integral_constant<unsigned long long, 1235>,
                    std::integral_constant<unsigned long long, 12345>,
                    std::integral_constant<unsigned long long, 123456>,
                    std::integral_constant<unsigned long long, 1234567>,
                    std::integral_constant<unsigned long long, 12345678>,
                    std::integral_constant<unsigned long long, 123456789>) {
    using Expr = typename RandomExpr<Seed::value, 5>::type;
    Expr expr{};
    auto dexpr = make_diff(expr, Variable<0>{});
    INFO("exrp=", expr.print());
    INFO("dexrp=", dexpr.print());
    Real x0 = Real(1.0);
    Real numeric = central_difference(expr, x0);
    CHECK(dexpr(x0) == doctest::Approx(numeric).epsilon(1e-6));
}
