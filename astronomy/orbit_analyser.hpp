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
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Difference;
using quantities::Length;
using quantities::MeasurementResult;
using quantities::Time;

// A characterization of the deviation of a value from some nominal value.
// The |half_width| is a 95th percentile, that is, this describes a value which
// is within [nominal_value - half_width, nominal_value + half_width] 95% of the
// time.
template<typename T>
struct Variability {
  T nominal_value;
  Difference<T> half_width;
};

template<typename T>
Variability<T> SequenceVariability(std::vector<T> const& values,
                                   T const& nominal_value);

template<typename Frame>
class OrbitAnalyser {
 public:
  OrbitAnalyser(not_null<Ephemeris<Frame>*> ephemeris,
                not_null<RotatingBody<Frame> const*> primary,
                not_null<RotatingBody<Frame> const*> sun,
                Instant initial_time,
                DegreesOfFreedom<Frame> initial_state,
                std::string name);

  OrbitAnalyser(not_null<Ephemeris<Frame>*> ephemeris,
                not_null<RotatingBody<Frame> const*> primary,
                not_null<RotatingBody<Frame> const*> sun,
                DiscreteTrajectory<Frame> const& trajectory,
                std::string name);

  void Analyse();

  void RecomputeProperties();

 private:
  not_null<Ephemeris<Frame>*> ephemeris_;
  not_null<RotatingBody<Frame> const*> primary_;
  not_null<RotatingBody<Frame> const*> sun_;
  DiscreteTrajectory<Frame> trajectory_;
  std::string name_;

  Instant reference_perihelion_time_;
  MeasurementResult<Time> tropical_year_;
  MeasurementResult<Angle> longitude_of_perihelion_;

  MeasurementResult<Time> nodal_period_;
  MeasurementResult<Time> anomalistic_period_;
  MeasurementResult<Time> sidereal_period_;
  MeasurementResult<AngularFrequency> nodal_precession_;
  MeasurementResult<AngularFrequency> apsidal_precession_;
  MeasurementResult<Length> periapsis_distance_;
  MeasurementResult<Length> apoapsis_distance_;
  MeasurementResult<double> eccentricity_;
  MeasurementResult<Angle> inclination_;
};

}  // namespace internal_orbit_analyser

using internal_orbit_analyser::OrbitAnalyser;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/orbit_analyser_body.hpp"
