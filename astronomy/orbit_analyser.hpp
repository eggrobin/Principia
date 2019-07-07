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

// All references in this class are to Capderou (2012), Satellites : de Kepler
// au GPS; the notation follows that of the book.
//
// The ground track of a satellite recurs when, in a reference frame fixing the
// plane of the orbit and the equator of the Earth (such a reference frame is
// possible under the assumption that the orbital plane precesses around the
// equator), the following two conditions occur simultaneously:
// — the satellite has gone through an integer number of revolutions Nᴛₒ;
// — the primary has gone through an integer number of revolutions Cᴛₒ.
// Such a recurrence is characterized by the value Nᴛₒ / Cᴛₒ.
// This class represents that rational number; while it exposes many quantities
// that are relevant to mission design, they are all derived from it.
//
// In this class, the word “day” is used to denote a revolution of the primary
// in the reference frame fixing the orbital plane.  Note that this day depends
// on the nodal precession of the satellite: it is not any of the usual
// definitions of the day:
// — for a sun-synchronous orbit (see sections 7.5.4, 7.5.5), this day is the
//   mean solar day;
// — for an orbit with no nodal precession (, i.e., a strictly polar orbit, see
//   section 7.1.3), this day is the stellar day.
// These days are counted negatively for a body with retrograde rotation.
class OrbitRecurrence final {
 public:
  // The following conditions must hold:
  //   Cᴛₒ ≠ 0;
  //   sign Cᴛₒ = sign νₒ if νₒ ≠ 0;
  //   |Dᴛₒ / Cᴛₒ| ≤ 1/2;
  //   gcd(Dᴛₒ, Cᴛₒ) = 1.
  OrbitRecurrence(int νₒ, int Dᴛₒ, int Cᴛₒ);

  template<typename Frame>
  static OrbitRecurrence ClosestRecurrence(
      Time const& nodal_period,
      AngularFrequency const& nodal_precession,
      RotatingBody<Frame> const& primary,
      int max_abs_Cᴛₒ);

  // The Capderou recurrence triple [νₒ; Dᴛₒ; Cᴛₒ], see section 11.1.3.
  // This triple expresses the rational Nᴛₒ / Cᴛₒ as an integer and fractional
  // part, Nᴛₒ / Cᴛₒ = νₒ + Dᴛₒ/ Cᴛₒ.

  // While he does mention orbits around bodies whose rotation is retrograde
  // (Venus and Triton), Capderou does not cover ground track recurrence in that
  // case.  We derive the sign from equations (11.4) and (11.5): Nᴛₒ is always
  // positive, and Cᴛₒ is negative for a primary with retrograde rotation.

  // νₒ is the number of orbits per day, rounded to the nearest integer.
  // Note that νₒ < 0 if the rotation of the primary is retrograde.
  int νₒ() const;
  int Dᴛₒ() const;
  // Cᴛₒ is the length of the cycle in days.
  // Note that Cᴛₒ < 0 if the rotation of the primary is retrograde.
  int Cᴛₒ() const;

  // The number of nodal periods per cycle Nᴛₒ, see section 11.1.2.
  int number_of_revolutions() const;

  // The equatorial shift Δλᴇ, see section 8.3.2.
  // This is counted eastward; it is negative (westward shift) except when
  // orbiting bodies whose rotation is retrograde.
  Angle equatorial_shift() const;

  // The width δʀ of the base interval, see 11.5.2.  This is always positive (it
  // represents the spacing between consecutive ground tracks on the equator).
  // Note that δʀ = |Δλᴇ|.
  Angle base_interval() const;
  // The grid interval δ, see section 11.5.2.
  // This is always positive (it represents the closest spacing between ground
  // tracks on the equator).
  Angle grid_interval() const;

  // The subcycle Eᴛₒ*, see section 11.5.3.
  // After Eᴛₒ* days, the ground track passes δ away from the origin.
  // Note that Eᴛₒ* < 0 if the rotation of the primary is retrograde.
  int subcycle() const;

 private:
  int νₒ_;
  int Dᴛₒ_;
  int Cᴛₒ_;
  int subcycle_;
};

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
