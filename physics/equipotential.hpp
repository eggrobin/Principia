#pragma once

#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_equipotential {

using base::not_null;
using geometry::Bivector;
using geometry::InfiniteFuture;
using geometry::Instant;
using geometry::Position;
using integrators::AdaptiveStepSizeIntegrator;
using integrators::ExplicitFirstOrderOrdinaryDifferentialEquation;
using physics::Ephemeris;
using quantities::Acceleration;
using quantities::Difference;
using quantities::Length;
using quantities::Speed;
using quantities::si::Metre;
using quantities::si::Second;

// TODO(phl): Similar class in Ephemeris.  Move to a common place.
template<typename ODE>
class ODEAdaptiveStepParameters final {
  public:
  ODEAdaptiveStepParameters(AdaptiveStepSizeIntegrator<ODE> const& integrator,
                            std::int64_t max_steps,
                            Length const& length_integration_tolerance);

  AdaptiveStepSizeIntegrator<ODE> const& integrator() const;
  std::int64_t max_steps() const;
  Length length_integration_tolerance() const;

  private:
  // This will refer to a static object returned by a factory.
  not_null<AdaptiveStepSizeIntegrator<ODE> const*> integrator_;
  std::int64_t max_steps_;
  Length length_integration_tolerance_;
};

template<typename Frame>
class Equipotential {
  static_assert(Frame::is_inertial);
 public:

  // The first state variable is a point of an equipotential.  The second state
  // variable is the (scale-free) intensity of the braking.
  using ODE =
      ExplicitFirstOrderOrdinaryDifferentialEquation<Position<Frame>, double>;
  using AdaptiveParameters = ODEAdaptiveStepParameters<ODE>;

  Equipotential(AdaptiveParameters const& adaptive_parameters,
                Ephemeris<Frame> const& ephemeris);

  std::vector<Position<Frame>> ComputeLine(
      Bivector<double, Frame> const& plane,
      Position<Frame> const& position,
      Instant const& t);

 private:
  // TODO(phl): Avoid Instant.
  using IndependentVariable = Instant;
  static constexpr IndependentVariable const s_initial_;
  static constexpr IndependentVariable const s_final_ = InfiniteFuture;
  static constexpr Difference<IndependentVariable> const s_initial_step_ =
      1 * Second;
  static constexpr Acceleration const characteristic_acceleration_ =
      1 * Metre / Second / Second;
  static constexpr Speed const characteristic_speed_ = 1 * Metre / Second;

  using State = typename ODE::State;
  using StateVariation = typename ODE::StateVariation;
  using SystemState = typename ODE::SystemState;
  using SystemStateError = typename ODE::SystemStateError;

  absl::Status RightHandSide(Bivector<double, Frame> const& plane,
                             Position<Frame> const& position,
                             Instant const& t,
                             IndependentVariable const& s,
                             State const& state,
                             StateVariation& state_variation);

  double ToleranceToErrorRatio(
      Difference<IndependentVariable> const& current_s_step,
      SystemStateError const& error);

  AdaptiveParameters const& adaptive_parameters_;
  not_null<Ephemeris<Frame> const*> const ephemeris_;
};

}  // namespace internal_equipotential

using internal_equipotential::Equipotential;

}  // namespace physics
}  // namespace principia

#include "physics/equipotential_body.hpp"
