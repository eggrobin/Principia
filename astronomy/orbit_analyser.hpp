#pragma once

#include "base/not_null.hpp"
#include "physics/ephemeris.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbit_analyser {

using base::not_null;
using geometry::Instant;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::RotatingBody;
using quantities::AngularFrequency;
using quantities::Time;

template<typename Frame>
class OrbitAnalyser {
 public:
  OrbitAnalyser(not_null<Ephemeris<Frame>*> ephemeris,
                not_null<RotatingBody<Frame> const*> primary,
                Instant initial_time,
                DegreesOfFreedom<Frame> initial_state);

  void Analyse();

 private:
  void RecomputeProperties();

  not_null<Ephemeris<Frame>*> ephemeris_;
  not_null<RotatingBody<Frame> const*> primary_;
  DiscreteTrajectory<Frame> trajectory_;

  Time nodal_period_;
  Time anomalistic_period_;
  Time sidereal_period_;
  AngularFrequency nodal_precession_;
  AngularFrequency apsidal_precession_;
};

}  // namespace internal_orbit_analyser

using internal_orbit_analyser::OrbitAnalyser;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/orbit_analyser_body.hpp"
