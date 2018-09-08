
#pragma once

#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <optional>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_generalized_runge_kutta_nyström_integrator {  // NOLINT(whitespace/line_length)

using base::make_not_null_unique;
using geometry::Sign;
using numerics::DoublePrecision;
using quantities::DebugString;
using quantities::Difference;
using quantities::Quotient;

template<typename Method, typename Position>
EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, Position>::
EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator() {
  // The first node is always 0 in an explicit method.
  CHECK_EQ(0.0, c_[0]);
  if (first_same_as_last) {
    // Check that the conditions for the FSAL property are satisfied, see for
    // instance Dormand, El-Mikkawy and Prince (1986),
    // Families of Runge-Kutta-Nyström formulae, equation 3.1.
    CHECK_EQ(1.0, c_[stages_ - 1]);
    CHECK_EQ(0.0, b̂_[stages_ - 1]);
    for (int j = 0; j < stages_ - 1; ++j) {
      CHECK_EQ(b̂_[j], a_[stages_ - 1][j]);
    }
  }
}

template<typename Method, typename Position>
Status EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
    Method,
    Position>::Instance::Solve(Instant const& t_final) {
}

template<typename Method, typename Position>
EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, Position> const&
EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, Position>::
Instance::integrator() const {
  return integrator_;
}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    ExplicitSecondOrderOrdinaryDifferentialEquation<Position>>::Instance>>
EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, Position>::
Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename Position>
void EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, Position>::
Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  AdaptiveStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  auto const& rkng_extension = serialization::
      EmbeddedExplicitGeneralizedRungeKuttaNystromIntegratorInstance::extension;
  auto* const extension =
      message
          ->MutableExtension(
              serialization::AdaptiveStepSizeIntegratorInstance::extension)
          ->MutableExtension(rkng_extension);
}

template<typename Method, typename Position>
EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, Position>::
Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio,
    Parameters const& parameters,
    EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator const& integrator)
    : AdaptiveStepSizeIntegrator<ODE>::Instance(
          problem, append_state, tolerance_to_error_ratio, parameters),
      integrator_(integrator) {}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    ExplicitSecondOrderOrdinaryDifferentialEquation<Position>>::Instance>>
EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, Position>::
NewInstance(IntegrationProblem<ODE> const& problem,
            AppendState const& append_state,
            ToleranceToErrorRatio const& tolerance_to_error_ratio,
            Parameters const& parameters) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(new Instance(problem,
                                                append_state,
                                                tolerance_to_error_ratio,
                                                parameters,
                                                *this));
}

template<typename Method, typename Position>
void EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, Position>::
WriteToMessage(not_null<serialization::AdaptiveStepSizeIntegrator*> message)
    const {
  message->set_kind(Method::kind);
}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    ExplicitSecondOrderOrdinaryDifferentialEquation<Position>>::Instance>>
EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method, Position>::
ReadFromMessage(
    serialization::AdaptiveStepSizeIntegratorInstance const& message,
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio,
    Parameters const& parameters) const {
  CHECK(message.HasExtension(
      serialization::
          EmbeddedExplicitGeneralizedRungeKuttaNystromIntegratorInstance::
              extension))
      << message.DebugString();

  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<typename AdaptiveStepSizeIntegrator<ODE>::Instance>(
      new Instance(problem,
                   append_state,
                   tolerance_to_error_ratio,
                   parameters,
                   *this));
}

}  // namespace internal_embedded_explicit_generalized_runge_kutta_nyström_integrator  // NOLINT

template<typename Method, typename Position>
internal_embedded_explicit_generalized_runge_kutta_nyström_integrator::
    EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method,
                                                           Position> const&
        EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator() {
  static_assert(
      std::is_base_of<methods::EmbeddedExplicitGeneralizedRungeKuttaNyström,
                      Method>::value,
      "Method must be derived from "
      "EmbeddedExplicitGeneralizedRungeKuttaNyström");
  static internal_embedded_explicit_generalized_runge_kutta_nyström_integrator::
      EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<Method,
                                                             Position> const
          integrator;
  return integrator;
}

}  // namespace integrators
}  // namespace principia
