#pragma once

#include "base/not_null.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/uncertainty.hpp"

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
using quantities::Length;
using quantities::MeasurementResult;
using quantities::Time;

template<typename Frame>
class OrbitAnalyser {
 public:
  OrbitAnalyser(not_null<Ephemeris<Frame>*> ephemeris,
                not_null<RotatingBody<Frame> const*> primary,
                Instant initial_time,
                DegreesOfFreedom<Frame> initial_state,
                std::string name);

  OrbitAnalyser(not_null<Ephemeris<Frame>*> ephemeris,
                not_null<RotatingBody<Frame> const*> primary,
                DiscreteTrajectory<Frame> const& trajectory,
                std::string name);

  void Analyse();

  void RecomputeProperties();

 private:
  not_null<Ephemeris<Frame>*> ephemeris_;
  not_null<RotatingBody<Frame> const*> primary_;
  DiscreteTrajectory<Frame> trajectory_;
  std::string name_;

  MeasurementResult<Time> nodal_period_;
  MeasurementResult<Time> anomalistic_period_;
  MeasurementResult<Time> sidereal_period_;
  MeasurementResult<AngularFrequency> nodal_precession_;
  MeasurementResult<AngularFrequency> apsidal_precession_;
  MeasurementResult<Length> periapsis_distance_;
  MeasurementResult<Length> apoapsis_distance_;
};

}  // namespace internal_orbit_analyser

using internal_orbit_analyser::OrbitAnalyser;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/orbit_analyser_body.hpp"
