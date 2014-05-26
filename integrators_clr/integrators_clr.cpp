#include "integrators_clr/integrators_clr.hpp"

#include <algorithm>
#include <vector>

#include "geometry/sign.hpp"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "testing_utilities/numerics.hpp"

using principia::testing_utilities::AbsoluteError;

namespace principia {
namespace integrators {

namespace {

inline void compute_harmonic_oscillator_force(double const t,
                                              std::vector<double> const& q,
                                              std::vector<double>* result) {
  (*result)[0] = -q[0];
}

inline void compute_harmonice_oscillator_velocity(std::vector<double> const& p,
                                                  std::vector<double>* result) {
  (*result)[0] = p[0];
}

}  // namespace


void SPRKTests::Run(System::Double tmax) {
  SPRKIntegrator integrator_;
  SPRKIntegrator::Parameters parameters_;
  SPRKIntegrator::Solution solution_;
  parameters_.q0 = {1.0};
  parameters_.p0 = {0.0};
  parameters_.t0 = 0.0;
#ifdef _DEBUG
  parameters_.tmax = 100.0;
#else
  parameters_.tmax = tmax;
#endif
  parameters_.Δt = 1.0E-4;
  parameters_.coefficients = integrator_.Order5Optimal();
  parameters_.sampling_period = 0;
  integrator_.Solve(&compute_harmonic_oscillator_force,
                         &compute_harmonice_oscillator_velocity,
                         parameters_,
                         &solution_);
  double q_error = 0;
  double p_error = 0;
  for (size_t i = 0; i < solution_.time.quantities.size(); ++i) {
    q_error = std::max(q_error,
                       std::abs(solution_.position[0].quantities[i] -
                                std::cos(solution_.time.quantities[i])));
    p_error = std::max(p_error,
                       std::abs(solution_.momentum[0].quantities[i] +
                                std::sin(solution_.time.quantities[i])));
  }
}

}  // namespace integrators
}  // namespace principia
