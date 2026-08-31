#include <auto_diff.hpp>
#include <auto_diff_base.hpp>
#include <auto_diff_squise.hpp>
#include <auto_diff_subst.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <fstream>
#include <ham.hpp>
#include <iostream>
#include <newton.hpp>
#include <numbers>
#include <parallel.hpp>
#include <precision.hpp>
#include <runge_kutta.hpp>
using namespace boost::numeric::ublas;

const Real e = 0.12;
const Real a = 1;
const Real n = pow(a, -1.5);
vector<Real> cache = vector<Real>(10000);
inline Real distance_sqr(Real t) {
  auto M = n * t;
  auto E = M;

  for (auto i = 0; i < 50; i++)
    E = E - (E - M - e * sin(E)) / (1 - e * cos(E));

  auto r = a * (1 - e * cos(E));
  return r * r;
}

inline Real distance_sqr_ch(Real t) {
  if (t >= 2 * std::numbers::pi)
    t = t - 2 * std::numbers::pi;
  auto M = n * t / (2 * std::numbers::pi) * (cache.size() - 1);
  size_t r = (size_t)floor(M);
  auto a = M - r;

  return cache[r] * (1 - a) + cache[r + 1] * a;
}

int main() {
  Parallel::for_i(
      [&](size_t i) -> void {
        cache[i] = distance_sqr(2 * std::numbers::pi / (cache.size() - 1) * i);
      },
      0, cache.size());

  auto q = Variable<0>{};
  auto p = Variable<1>{};
  auto t = Variable<2>{};
  auto r = apply_function<distance_sqr_ch>(t);
  auto z = sqrt(r + q * q);
  auto ham = p * p * Constant<0.5>{} - Constant<1.0>{} / z;

  auto rhs = Hamiltonian::compile_system<1>(ham);
  int N = 128;
  int NSteps = 1000;
  auto inits = std::vector<vector<Real>>(N * N);
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      inits[i + j * N] = vector<Real>(2);
      inits[i + j * N][0] = -2.5 + 5.0 / (N - 1) * i;
      inits[i + j * N][1] = -1 + 2.0 / (N - 1) * j;
    }
  }

  std::ofstream f("out.dat");

  for (auto i = 0; i < NSteps; i++) {
    Parallel::for_i(
        [&](size_t i) -> void {
          auto diff_solver = RungeKutta::BRK(2, 3);
          RungeKutta::Integrate(2, inits[i], 0, 2 * std::numbers::pi, 0.1, rhs,
                                diff_solver);
        },
        0, N *N);
    std::cerr << std::format("\r{:.2f}%", (Real)i / (NSteps - 1) * 100);
    for (size_t j = 0; j < inits.size(); j++) {
      for (auto &v : inits[j]) {
        f << v << " ";
      }
      f << '\n';
    }
  }

  return 0;
}
