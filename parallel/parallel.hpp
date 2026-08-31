#pragma once
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <omp.h>
#include <utility>
namespace Parallel {
template <typename U, typename T>
concept IndexableContainer = requires(U u, size_t i) {
  { u[i] } -> std::same_as<T &>;
  { u[0] } -> std::same_as<T &>;
  { u.size() } -> std::same_as<size_t>;
};

template <typename Op, typename T>
concept BinaryOperation = requires(Op &&op, const T &a, const T &b) {
  { op(a, b) } -> std::same_as<T>;
};
namespace detail {
template <typename T, typename Op, typename U>
  requires BinaryOperation<Op, T> && IndexableContainer<U, T>
T reduce_seq(Op &&binary_op, const U &collection, size_t start, size_t end,
             T identity, int depth, int max_depth) {
  if (depth > max_depth || (end - start) < 128) {
    T result = identity;
    for (size_t i = start; i < end; i++) {
      result = binary_op(result, collection[i]);
    }
    return result;
  }

  T left, right;
  size_t mid = start + (end - start) / 2;
#pragma omp task shared(left)
  {
    left = reduce_seq(std::forward<Op>(binary_op), collection, start, mid,
                      identity, depth + 1, max_depth);
  }

  right = reduce_seq(std::forward<Op>(binary_op), collection, mid, end,
                     identity, depth + 1, max_depth);

#pragma omp taskwait

  return binary_op(left, right);
}
} // namespace detail

template <typename T, typename Op, typename U>
  requires BinaryOperation<Op, T> && IndexableContainer<U, T>
T reduce(Op &&binary_op, const U &collection, T identity) {
  T result = identity;
#pragma omp parallel
  {
#pragma omp single nowait
    {
      int n = omp_get_max_threads();
      int max_depth = ceil(log2(n));
      result = detail::reduce_seq(std::forward<Op>(binary_op), collection, 0,
                                  collection.size(), identity, 1, max_depth);
    }
  }
  return result;
}

template <typename Fn>
  requires std::invocable<Fn &, size_t>
void for_i(Fn &&f, int start, int end) {
  int N = std::abs(end - start);
  int max_threads = std::max(std::min(N, omp_get_max_threads()), 1);

#pragma omp parallel for num_threads(max_threads) schedule(dynamic, 10)
  for (size_t i = start; i < end; i++) {
    std::invoke(std::forward<Fn>(f), i);
  }
}

template <typename T, typename Fn, typename U>
  requires std::invocable<Fn &, const T &> && IndexableContainer<U, T>
void for_each(Fn &&f, U &collection) {
  int N = collection.size();
  int max_threads = std::max(std::min(N, omp_get_max_threads()), 1);

#pragma omp parallel for num_threads(max_threads) schedule(dynamic, 10)
  for (size_t i = 0; i < N; i++) {
    std::invoke(std::forward<Fn>(f), collection[i]);
  }
}
} // namespace Parallel
