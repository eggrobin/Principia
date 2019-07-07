﻿#pragma once

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
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Difference;
using quantities::Length;
using quantities::Time;

// A characterization of the deviation of a value from some nominal value.
// The result is a 95th percentile, that is, with
//   half_width = Variability(values, nominal_value),
// 95% of the |values| are within
// [nominal_value - half_width, nominal_value + half_width].
template<typename T>
Difference<T> Variability(std::vector<T> const& values, T const& nominal_value);

template<typename Frame>
class OrbitAnalyser final {
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
  Time tropical_year_;
  Angle longitude_of_perihelion_;

  Time nodal_period_;
  Time anomalistic_period_;
  Time sidereal_period_;
  AngularFrequency nodal_precession_;
  AngularFrequency apsidal_precession_;
  Length periapsis_distance_;
  Length apoapsis_distance_;
  double eccentricity_;
  Angle inclination_;
};

}  // namespace internal_orbit_analyser

using internal_orbit_analyser::OrbitAnalyser;

}  // namespace astronomy
}  // namespace principia

#include "astronomy/orbit_analyser_body.hpp"
