#pragma once

#include "integrators/integrators.hpp"
#include "integrators/ordinary_differential_equations.hpp"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_runge_kutta_instance {

using base::Status;
using geometry::Instant;

template<typename ODE, typename Method>
class EmbeddedExplicitRungeKuttaInstance
    : public typename AdaptiveStepSizeIntegrator<ODE>::Instance {
  using typename Integrator<ODE>::AppendState;
  using typename AdaptiveStepSizeIntegrator<ODE>::Parameters;
  using typename AdaptiveStepSizeIntegrator<ODE>::ToleranceToErrorRatio;
 public:
  Status Solve(Instant const& t_final) override;
};

}  // namespace internal_embedded_explicit_runge_kutta_instance
}  // namespace integrators
}  // namespace principia
