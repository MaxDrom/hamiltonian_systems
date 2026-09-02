#include <auto_diff.hpp>
#include <auto_diff_base.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
#include <ham.hpp>
#include <numbers>
#include <parallel.hpp>
#include <vector>
#include <yoshida/leapfrog.hpp>
#include <yoshida/leapfrog_extension.hpp>
#include <yoshida/yoshida.hpp>
using namespace boost::numeric::ublas;
using namespace AutoDiff;
int main() {
    auto q = Variable<0>{};
    auto p = Variable<1>{};
    constexpr double omega = 2 * std::numbers::pi;
    auto ham = p * p * Constant<0.5>{} - cos(Constant<omega>{} * q);

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
            auto leapfrog = Yoshida::make_leapfrog_for_sep<1>(ham);
            auto integrator = Yoshida::make_yoshida_from_leapfrog<4>(leapfrog);
            RungeKutta::Integrate(inits[i], 0, 2 * std::numbers::pi, 0.001,
                                  integrator,
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
