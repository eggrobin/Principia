﻿
#pragma once

#include <experimental/optional>
#include <memory>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/forkable.hpp"
#include "physics/massive_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {

namespace geometry {
template <typename FrameTag, FrameTag frame_tag, bool frame_is_inertial> class Frame;
template <typename Scalar, typename Frame, int rank> class Multivector;
}  // namespace geometry

using geometry::Instant;
using geometry::Vector;
using physics::DiscreteTrajectory;
using physics::DynamicFrame;
using physics::Ephemeris;
using physics::Frenet;
using quantities::Force;
using quantities::Mass;
using quantities::SpecificImpulse;
using quantities::Speed;
using quantities::Time;
using quantities::Variation;

namespace ksp_plugin {

// This class represents a constant-thrust inertial burn.  |InertialFrame| is
// an underlying inertial reference frame, |Frame| is the reference frame used
// to compute the Frenet frame.  |Frame| is defined by the parameter |frame|
// given to the constructor.  The |direction| is given in the Frenet frame of
// the trajectory at the beginning of the burn.
template<typename InertialFrame, typename Frame>
class Manœuvre {
 public:
  Manœuvre(Force const& thrust,
           Mass const& initial_mass,
           SpecificImpulse const& specific_impulse,
           Vector<double, Frenet<Frame>> const& direction,
           not_null<std::unique_ptr<DynamicFrame<InertialFrame, Frame> const>>
               frame);
  Manœuvre(Manœuvre&&) = default;
  Manœuvre& operator=(Manœuvre&&) = default;
  ~Manœuvre() = default;

  Force const& thrust() const;
  Mass const& initial_mass() const;
  // Specific impulse by mass, because specific impulse by weight is insane.
  // This is defined as the ratio of thrust to mass flow.
  // If the burn is done with a single engine (in a vacuum), this will be its
  // exhaust velocity.  For several engines, this is the total thrust divided
  // by the sum of the individual mass flows (where each mass flow is the
  // individual thrust divided by the exhaust velocity).
  SpecificImpulse const& specific_impulse() const;
  Vector<double, Frenet<Frame>> const& direction() const;
  not_null<DynamicFrame<InertialFrame, Frame> const*> frame() const;

  // Equivalent characterizations of intensity.  Only one of the mutators may be
  // called, and only once.
  Time duration() const;
  void set_duration(Time const& duration);
  Speed Δv() const;
  void set_Δv(Speed const& Δv);

  // Equivalent characterizations of timing.  Only one of the mutators may be
  // called, and only once.
  Instant initial_time() const;
  void set_initial_time(Instant const& initial_time);
  // Intensity and timing must have been set.
  Instant time_of_half_Δv() const;
  // |Δv| or |duration| must have been set.
  void set_time_of_half_Δv(Instant const& time_of_half_Δv);

  // Derived quantities.
  Variation<Mass> mass_flow() const;
  // Intensity must have been set.
  Mass final_mass() const;
  // Intensity must have been set.
  Time time_to_half_Δv() const;
  // Intensity and timing must have been set.
  Instant final_time() const;

  // Intensity and timing must have been set.
  // Returns true if and only if [initial_time, final_time] ⊆ ]begin, end[.
  bool FitsBetween(Instant const& begin, Instant const& end) const;

  // Intensity and timing must have been set.  The result is valid until
  // |*this| is destroyed.  |coasting_trajectory| must have a point at
  // |initial_time()|.
  typename Ephemeris<InertialFrame>::IntrinsicAcceleration acceleration(
      DiscreteTrajectory<InertialFrame> const& coasting_trajectory) const;

 private:
  Force const thrust_;
  Mass const initial_mass_;
  SpecificImpulse const specific_impulse_;
  Vector<double, Frenet<Frame>> const direction_;
  std::experimental::optional<Time> duration_;
  std::experimental::optional<Instant> initial_time_;
  not_null<std::unique_ptr<DynamicFrame<InertialFrame, Frame> const>> frame_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/manœuvre_body.hpp"
