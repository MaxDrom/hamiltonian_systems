#include "kinetic_integrator.hpp"
#include "triangle.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
#include <numbers>
#include <runge_kutta/odesolver.hpp>
#include <yoshida/leapfrog_extension.hpp>
#include <yoshida/yoshida.hpp>
int main() {
    auto V = get_potential<1.0, 1.0, 1.0>();
    auto leapfrog = Yoshida::make_leapfrog_for_spq<3>(V, kinetic_solver);
    auto integrator = Yoshida::make_yoshida_from_leapfrog<6>(leapfrog);

    using namespace boost::numeric::ublas;
    vector<Real> inits(6);

    inits(0) = 2;
    inits(1) = 0;
    inits(2) = 0;

    inits(3) = 0;
    inits(4) = 2.4048963567502017 / 8;
    inits(5) = 3.690826688109932 / 8;
    auto f = std::ofstream("out.dat");
    RungeKutta::Integrate(inits, 0, 4 * std::numbers::pi, 0.001, integrator,
                          [&f](Real t, const vector<Real> &x) -> void {
                              auto r = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
                              r = mp_sqrt(r);
                              for (auto &v : x)
                                  f << v / r << " ";
                              f << "\n";
                          });
}
