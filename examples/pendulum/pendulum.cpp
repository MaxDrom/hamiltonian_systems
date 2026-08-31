#include <auto_diff.hpp>
#include <auto_diff_base.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
#include <ham.hpp>
#include <numbers>
#include <parallel.hpp>
#include <runge_kutta.hpp>
#include <vector>
using namespace boost::numeric::ublas;

int main() {
    auto q = Variable<0>{};
    auto p = Variable<1>{};
    constexpr double omega = 2 * std::numbers::pi;
    auto ham = p * p * Constant<0.5>{} - cos(Constant<omega>{} * q);
    auto rhs = Hamiltonian::compile_system<1>(ham);

    int N = 20;
    auto inits = std::vector<vector<Real>>(N);

    for (size_t i = 0; i < N / 2; i++) {
        inits[i] = vector<Real>(2);
        inits[i][0] = 0.5 * i / (N / 2.0 - 1);
        inits[i][1] = 0;
    }

    for (size_t i = N / 2; i < N; i++) {
        inits[i] = vector<Real>(2);
        inits[i][1] = -2.0 + 4.0 * (i - N / 2.0) / (N / 2.0 - 1);
        inits[i][0] = inits[i][1] < 0 ? 0.5 : -0.5;
    }

    auto results = std::vector<std::vector<vector<Real>>>(N);
    Parallel::for_i(
        [&](size_t i) -> void {
            auto integrator = RungeKutta::BRK(2, 3);
            RungeKutta::Integrate(2, inits[i], 0, 2 * std::numbers::pi, 0.001,
                                  rhs, integrator,
                                  [&](Real t, const vector<Real> &x) -> void {
                                      results[i].push_back(x);
                                  });
        },
        0, N);

    std::ofstream f("traj.dat");

    for (size_t i = 0; i < results[0].size(); i++) {
        for (size_t j = 0; j < N; j++) {
            for (auto &v : results[j][i])
                f << v << " ";
        }
        f << '\n';
    }
}
