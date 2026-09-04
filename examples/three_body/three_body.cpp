#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
#include <ham.hpp>
#include <ham_maker.hpp>
#include <numbers>
#include <precision.hpp>
#include <runge_kutta/brk.hpp>
using namespace boost::numeric::ublas;
using namespace AutoDiff;
int main() {
    auto body1 = Body<1.0>{};
    auto body2 = Body<1.0>{};
    auto body3 = Body<1.0>{};
    auto problem = make_planar_nbody(body1, body2, body3);
    auto system = problem.compiled_system();
    auto integrator = RungeKutta::make_brk<problem.NCoords>(2, 3, system);
    auto inits = vector<Real>(problem.NCoords);
    inits[0] = 0.97000436;  // x1
    inits[1] = -0.24308753; // y1
    inits[2] = 0.0;         // x2
    inits[3] = 0.0;         // y2
    inits[4] = -0.97000436; // x3
    inits[5] = 0.24308753;  // y3

    inits[6] = 0.4662036850;  // px1
    inits[7] = 0.4323657300;  // py1
    inits[8] = -0.9324073700; // px2
    inits[9] = -0.8647314600; // py2
    inits[10] = 0.4662036850; // px3
    inits[11] = 0.4323657300; // py3
    std::ofstream f("traj.dat");

    RungeKutta::Integrate(inits, 0, 2.01 * std::numbers::pi / 3, 0.001,
                          integrator,
                          [&](Real t, const vector<Real> &x) -> void {
                              for (auto &v : x)
                                  f << v << " ";
                              f << '\n';
                          });
    return 0;
}
