#include "newton.hpp"
#include <Eigen/Dense>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <functional>
#include <iostream>
#include <omp.h>
#include <precision.hpp>
using namespace boost::numeric::ublas;
namespace Newton {
using namespace boost::multiprecision;
void jacobi(matrix<Real> &jac, const vector<Real> &x,
            std::function<void(const vector<Real> &, vector<Real> &)> f,
            const vector<Real> &f_x) {
  Real h = sqrt(std::numeric_limits<Real>::epsilon());

  int N = x.size();
  int max_threads = std::min(N, omp_get_max_threads());

#pragma omp parallel for num_threads(max_threads)
  for (size_t i = 0; i < N; i++) {
    auto x_copy = x;
    x_copy(i) += h;
    vector<Real> new_f(x.size());
    f(x_copy, new_f);
    for (size_t j = 0; j < x.size(); j++) {
      jac(j, i) = (new_f(j) - f_x(j));
    }
  }
  jac = jac / h;
}

int solve_linear(const matrix<Real> &A, const vector<Real> &b,
                 vector<Real> &x) {
  // Проверяем корректность размеров квадратной системы
  if (A.size1() != A.size2() || A.size1() != b.size()) {
    return -1; // Ошибка размерности
  }

  // Изменяем размер результирующего вектора x под размер b
  x.resize(b.size(), false);

  const int rows = static_cast<int>(A.size1());
  const int cols = static_cast<int>(A.size2());

  // Типы для маппинга матрицы (uBLAS по умолчанию использует RowMajor)
  using EigenMatrixMapConst =
      Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>>;
  // Типы для маппинга векторов
  using EigenVectorMapConst =
      Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>>;
  using EigenVectorMap = Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>>;

  // Отображаем память Boost в структуры Eigen без накладных расходов на
  // копирование
  EigenMatrixMapConst eigen_A(&A.data()[0], rows, cols);
  EigenVectorMapConst eigen_b(&b.data()[0], rows);
  EigenVectorMap eigen_x(&x.data()[0], rows);

  // Используем Полное Ортогональное Разложение (COD) для плохих матриц
  Eigen::CompleteOrthogonalDecomposition<
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
      cod(eigen_A);

  // Находим псевдорешение с минимальной нормой и пишем напрямую в память
  // вектора x
  eigen_x = cod.solve(eigen_b);

  // Метод COD успешно разрешает даже вырожденные ситуации,
  // поэтому возвращаем 0 как маркер успешного завершения расчета.
  return 0;
}

int solve_linear(const matrix<Real> &A, const matrix<Real> &b,
                 matrix<Real> &x) {
  matrix<Real> LU = A;
  permutation_matrix pm(A.size1());
  auto res = lu_factorize(LU, pm);
  if (res != 0)
    return res;
  x = b;
  lu_substitute(LU, pm, x);
  return 0;
}

int inverse(const matrix<Real> &A, matrix<Real> &inverse) {
  matrix<Real> id = identity_matrix<Real>(A.size1());
  return solve_linear(A, id, inverse);
}

int solve(vector<Real> &x,
          std::function<void(const vector<Real> &, vector<Real> &)> f,
          Real damping) {
  vector<Real> dx(x.size());
  matrix<Real> jac(x.size(), x.size());
  vector<Real> fx(x.size());
  Real h = 10 * sqrt(std::numeric_limits<Real>::epsilon());
  do {
    f(x, fx);
    jacobi(jac, x, f, fx);
    int res = solve_linear(jac, fx, dx);
    if (res != 0) {
      // std::cout<<jac<<'\n';
      if (norm_2(fx) <= std::numeric_limits<Real>::min())
        return 0;
      else
        return 1;
    }
    std::cout << norm_2(fx) << " " << fx << '\n';
    x -= damping * dx;

  } while (norm_2(dx) > 10 * std::max(h, norm_2(x) * h));

  return 0;
}

int solve(Real &x, std::function<Real(Real)> f) {
  Real h = 2 * sqrt(std::numeric_limits<Real>::epsilon());
  Real dx;
  Real fx;
  Real hh = sqrt(std::numeric_limits<Real>::epsilon());
  do {
    fx = f(x);
    auto df = (f(x + h) - fx) / h;
    if (abs(df) < std::numeric_limits<Real>::min() || df == 0)
      return 1;
    dx = -fx / df;
    x = x + dx;

  } while (abs(dx) > 10 * std::max(hh, abs(x) * hh));

  return 0;
}

int solve(Real &x, std::function<Real(Real)> f, std::function<Real(Real)> df) {

  Real dx;
  Real fx;
  do {
    fx = f(x);
    auto dfx = df(x);
    if (abs(dfx) < std::numeric_limits<Real>::min() || dfx == 0)
      return 1;
    dx = -fx / dfx;
    x = x + dx;
    // std::cout<<x<<" "<<dx<<'\n';
  } while (abs(fx) > std::numeric_limits<Real>::min() &&
           abs(dx) > std::numeric_limits<Real>::epsilon());

  return 0;
}

Real find_minimum(Real a, Real b, std::function<Real(Real)> f) {
  Real phi = (sqrt(Real(5)) - Real(1)) / 2;
  Real c, d, x;
  Real eps = std::numeric_limits<Real>::epsilon();
  do {
    c = b - (b - a) * phi;
    d = a + (b - a) * phi;
    if (f(c) < f(d))
      b = d;
    else
      a = c;

    x = (a + b) / 2;

  } while ((b - a) > std::max(eps, abs(x) * eps));

  return x;
}
} // namespace Newton
