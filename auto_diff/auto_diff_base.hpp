#pragma once
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <format>
#include <precision.hpp>
#include <string>

template <typename T>
concept IsReal = std::convertible_to<T, Real>;

template <typename Derived> struct Expression {
private:
  inline const Derived *pimpl() const {
    return static_cast<const Derived *>(this);
  }

public:
  template <IsReal... Args> inline Real operator()(Args &&...x) const {
    return pimpl()->apply_impl(std::forward<Args>(x)...);
  }

  inline std::string print() const { return pimpl()->print_impl(); }

  inline size_t priority() const { return pimpl()->priority_impl(); }
};

template <size_t Idx> struct Variable : public Expression<Variable<Idx>> {
  friend struct Expression<Variable<Idx>>;

private:
  template <IsReal... Args> inline Real apply_impl(Args &&...args) const {
    return std::get<Idx>(std::forward_as_tuple(std::forward<Args>(args)...));
  }
  inline std::string print_impl() const { return std::format("x_{}", Idx); }
  inline size_t priority_impl() const { return 100; }
};

template <typename T>
concept IsExpression = std::is_base_of_v<Expression<T>, T>;

template <typename T> struct is_variable_trait : std::false_type {};

template <std::size_t Index>
struct is_variable_trait<Variable<Index>> : std::true_type {};

template <typename T>
concept IsVariable = is_variable_trait<std::decay_t<T>>::value;

template <double Val> struct Constant : public Expression<Constant<Val>> {
  friend struct Expression<Constant<Val>>;

private:
  template <IsReal... Args> inline Real apply_impl(Args &&...args) const {
    return Val;
  }

  inline std::string print_impl() const { return std::format("{}", Val); }

  inline size_t priority_impl() const {
    if constexpr (Val < 0)
      return 0;
    return 100;
  }
};

template <int Val> struct Integer : public Expression<Integer<Val>> {
  friend struct Expression<Integer<Val>>;

private:
  template <IsReal... Args> inline Real apply_impl(Args &&...args) const {
    return Val;
  }

  inline std::string print_impl() const { return std::format("{}", Val); }

  inline size_t priority_impl() const {
    if constexpr (Val < 0)
      return 0;
    return 100;
  }
};

template <typename T> struct is_constant_trait : std::false_type {};

template <double Value>
struct is_constant_trait<Constant<Value>> : std::true_type {};

template <typename T>
concept IsConstant = is_constant_trait<std::decay_t<T>>::value;

template <typename T> struct is_integer_trait : std::false_type {};

template <int Value>
struct is_integer_trait<Integer<Value>> : std::true_type {};

template <typename T>
concept IsInteger = is_integer_trait<std::decay_t<T>>::value;

template <typename T>
concept IsNumber = IsInteger<T> || IsConstant<T>;

template <typename T>
concept IsNotNumber = !IsNumber<T>;