#pragma once
#include <auto_diff_base.hpp>
#include <precision.hpp>
#include <utility>

template <typename Derived> struct BinaryOpTag {
  private:
    inline const Derived *pimpl() const {
        return static_cast<const Derived *>(this);
    }

  public:
    inline Real apply(Real a, Real b) const {
        return pimpl()->apply_impl(a, b);
    }
    inline std::string get_char() const { return pimpl()->get_char_impl(); }

    inline size_t priority() const { return pimpl()->priority_impl(); }
};

template <typename Derived> struct FuntionTag {
  private:
    inline const Derived *pimpl() const {
        return static_cast<const Derived *>(this);
    }

  public:
    inline Real apply(Real &&a) const {
        return pimpl()->apply_impl(std::forward<Real>(a));
    }
    inline std::string get_string() const { return pimpl()->get_string_impl(); }
};

template <auto Func> constexpr std::string_view get_raw_function_name() {
#if defined(__clang__) || defined(__GNUC__)
    return __PRETTY_FUNCTION__;
#elif defined(_MSC_VER)
    return __FUNCSIG__;
#else
    return "unknown";
#endif
}

template <auto func>
struct FuntionTagFromFunc : FuntionTag<FuntionTagFromFunc<func>> {
    friend struct FuntionTag<FuntionTagFromFunc<func>>;

  private:
    inline Real apply_impl(Real &&a) const {
        return func(std::forward<Real>(a));
    }
    inline std::string get_string_impl() const {
        constexpr std::string_view raw = get_raw_function_name<func>();

        constexpr std::string_view target = "Func = ";
        constexpr auto pos = raw.find(target);

        if constexpr (pos == std::string_view::npos) {
            return std::string(raw);
        } else {

            constexpr auto start = pos + target.size();
            constexpr auto end = raw.find(']', start);
            std::string_view name = raw.substr(start, end - start);

            if (!name.empty() && name[0] == '&') {
                name.remove_prefix(1);
            }
            return std::string(name);
        }
    }
};

template <typename T>
concept IsFunctionTag = std::is_base_of_v<FuntionTag<T>, T>;

template <IsExpression A, IsFunctionTag FuncTag>
struct Function : Expression<Function<A, FuncTag>> {
    friend struct Expression<Function<A, FuncTag>>;

  private:
    [[no_unique_address]]
    A _value;
    [[no_unique_address]]
    FuncTag _func;

    template <IsReal... Args> inline Real apply_impl(Args &&...args) const {
        auto result = _value(std::forward<Args>(args)...);
        return _func.apply(std::forward<Real>(result));
    }
    inline std::string print_impl() const {
        return std::format("{}({})", _func.get_string(), _value.print());
    }

    inline size_t priority_impl() const { return 100; }
};

template <IsExpression Value>
using Sin =
    Function<Value, FuntionTagFromFunc<static_cast<Real (*)(Real &&)>(mp_sin)>>;

template <IsExpression Value> static inline Sin<Value> sin(const Value &v) {
    return {};
}

template <IsExpression Value>
using Cos =
    Function<Value, FuntionTagFromFunc<static_cast<Real (*)(Real &&)>(mp_cos)>>;

template <IsExpression Value> static inline Cos<Value> cos(const Value &v) {
    return {};
}

template <IsExpression Value>
using Sqrt =
    Function<Value,
             FuntionTagFromFunc<static_cast<Real (*)(Real &&)>(mp_sqrt)>>;

template <IsExpression Value> static inline Sqrt<Value> sqrt(const Value &v) {
    return {};
}

template <IsExpression Value>
using Exp =
    Function<Value, FuntionTagFromFunc<static_cast<Real (*)(Real &&)>(mp_exp)>>;

template <IsExpression Value> static inline Exp<Value> exp(const Value &v) {
    return {};
}

template <typename T>
concept IsBinaryOpTag = std::is_base_of_v<BinaryOpTag<T>, T>;

template <IsExpression Left, IsExpression Right, IsBinaryOpTag OpTag>
struct BinaryOp : public Expression<BinaryOp<Left, Right, OpTag>> {
    friend struct Expression<BinaryOp<Left, Right, OpTag>>;

  private:
    [[no_unique_address]]
    Left _left;
    [[no_unique_address]]
    Right _right;
    [[no_unique_address]]
    OpTag _op;
    template <IsReal... Args> inline Real apply_impl(Args &&...args) const {
        return _op.apply(_left(std::forward<Args>(args)...),
                         _right(std::forward<Args>(args)...));
    }
    inline std::string print_impl() const {
        auto left = _op.priority() >= _left.priority()
                        ? std::format("({})", _left.print())
                        : _left.print();
        auto right = _op.priority() >= _right.priority()
                         ? std::format("({})", _right.print())
                         : _right.print();
        return std::format("{}{}{}", left, _op.get_char(), right);
    }

    inline size_t priority_impl() const { return _op.priority(); }
};

template <IsExpression A, IsBinaryOpTag OpTag>
struct BinaryOp<A, A, OpTag> : Expression<BinaryOp<A, A, OpTag>> {
    friend struct Expression<BinaryOp<A, A, OpTag>>;

  private:
    [[no_unique_address]]
    A _left;
    [[no_unique_address]]
    OpTag _op;
    template <IsReal... Args> inline Real apply_impl(Args &&...args) const {
        Real value = _left(std::forward<Args>(args)...);
        return _op.apply(value, value);
    }
    inline std::string print_impl() const {
        auto left = _op.priority() > _left.priority()
                        ? std::format("({})", _left.print())
                        : _left.print();
        return std::format("{}{}{}", left, _op.get_char(), left);
    }
    inline size_t priority_impl() const { return _op.priority(); }
};

struct MulTag : public BinaryOpTag<MulTag> {
    friend struct BinaryOpTag<MulTag>;

  private:
    inline Real apply_impl(Real a, Real b) const { return a * b; }
    inline std::string get_char_impl() const { return "*"; }

    inline size_t priority_impl() const { return 2; }
};

struct AddTag : public BinaryOpTag<AddTag> {
    friend struct BinaryOpTag<AddTag>;

  private:
    inline Real apply_impl(Real a, Real b) const { return a + b; }

    inline std::string get_char_impl() const { return "+"; }
    inline size_t priority_impl() const { return 1; }
};

struct MinusTag : public BinaryOpTag<MinusTag> {
    friend struct BinaryOpTag<MinusTag>;

  private:
    inline Real apply_impl(Real a, Real b) const { return a - b; }

    inline std::string get_char_impl() const { return "-"; }

    inline size_t priority_impl() const { return 1; }
};

struct DivTag : public BinaryOpTag<DivTag> {
    friend struct BinaryOpTag<DivTag>;

  private:
    inline Real apply_impl(Real a, Real b) const { return a / b; }

    inline std::string get_char_impl() const { return "/"; }
    inline size_t priority_impl() const { return 2; }
};

template <IsExpression A, IsExpression B> using Mul = BinaryOp<A, B, MulTag>;

template <IsExpression A, IsExpression B>
constexpr static inline Mul<A, B> operator*(const Expression<A> &a,
                                            const Expression<B> &b) {
    return {};
};

template <IsExpression A, IsExpression B> using Div = BinaryOp<A, B, DivTag>;

template <IsExpression A, IsExpression B>
constexpr static inline Mul<A, Div<Constant<1.0>, B>>
operator/(const Expression<A> &a, const Expression<B> &b) {
    return {};
};

template <auto Func, IsExpression A>
constexpr static inline Function<A, FuntionTagFromFunc<Func>>
apply_function(A a) {
    return {};
}

template <IsExpression A, IsExpression B> using Add = BinaryOp<A, B, AddTag>;

template <IsExpression A, IsExpression B>
constexpr static inline Add<A, B> operator+(const Expression<A> &a,
                                            const Expression<B> &b) {
    return {};
};

template <IsExpression A, IsExpression B>
using Minus = BinaryOp<A, B, MinusTag>;

template <IsExpression A, IsExpression B>
constexpr static inline Add<A, Mul<Constant<-1.0>, B>>
operator-(const Expression<A> &a, const Expression<B> &b) {
    return {};
};

template <IsExpression A>
constexpr static inline Minus<Integer<0>, A> operator-(const Expression<A> &a) {
    return {};
};

template <template <typename...> class NodeWorker, IsExpression A,
          typename... Args>
struct ExpressionIterator {
    using type = typename NodeWorker<A, Args...>::type;
};

template <template <typename...> class NodeWorker, IsExpression Left,
          IsExpression Right, IsBinaryOpTag OpTag, typename... Args>
struct ExpressionIterator<NodeWorker, BinaryOp<Left, Right, OpTag>, Args...> {
    using _Left = typename ExpressionIterator<NodeWorker, Left, Args...>::type;
    using _Right =
        typename ExpressionIterator<NodeWorker, Right, Args...>::type;

    using type =
        typename NodeWorker<BinaryOp<_Left, _Right, OpTag>, Args...>::type;
};

template <template <typename...> class NodeWorker, IsExpression A,
          IsFunctionTag FuncTag, typename... Args>
struct ExpressionIterator<NodeWorker, Function<A, FuncTag>, Args...> {
    using SimplifiedV =
        typename ExpressionIterator<NodeWorker, A, Args...>::type;

    using type =
        typename NodeWorker<Function<SimplifiedV, FuncTag>, Args...>::type;
};

template <size_t N_iters, IsExpression Expr, IsExpression Pattern,
          IsExpression Initial>
constexpr decltype(auto) recursive_subst(Expr expr, Pattern pattern,
                                         Initial initial) {
    if constexpr (N_iters == 1)
        return simplify(subst(expr, pattern, initial));
    else {
        auto prev = recursive_subst<N_iters - 1>(expr, pattern, initial);
        return simplify(subst(expr, pattern, prev));
    }
}

template <size_t N_iters, IsExpression Expr, IsExpression Init, size_t Idx>
constexpr decltype(auto) newton_expr(Expr expr, Variable<Idx> var, Init init) {
    auto df = make_diff(expr, var);
    auto expr1 = var - expr / df;
    return recursive_subst<N_iters>(expr1, var, init);
}
