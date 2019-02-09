﻿#include "astronomy/orbit_analyser.hpp"

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"


namespace principia {

using astronomy::ICRS;
using astronomy::J2000;
using base::dynamic_cast_not_null;
using base::OFStream;
using geometry::Displacement;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MasslessBody;
using physics::OblateBody;
using physics::RelativeDegreesOfFreedom;
using physics::SolarSystem;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::ArcSin;
using quantities::Cbrt;
using quantities::Cos;
using quantities::Length;
using quantities::Pow;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::AstronomicalUnit;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;

namespace astronomy {

class OrbitAnalyzerTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    ephemeris_ = solar_system_2000_.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<ICRS>>(),
            /*step=*/10 * Minute));
  }

  static SolarSystem<ICRS> solar_system_2000_;
  static std::unique_ptr<Ephemeris<ICRS>> ephemeris_;
};

SolarSystem<ICRS> OrbitAnalyzerTest::solar_system_2000_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2451545_000000000.proto.txt");
std::unique_ptr<Ephemeris<ICRS>> OrbitAnalyzerTest::ephemeris_;

#if !defined(_DEBUG)

TEST_F(OrbitAnalyzerTest, DISABLED_Молния) {
  auto const earth_body = dynamic_cast_not_null<OblateBody<ICRS> const*>(
      solar_system_2000_.massive_body(*ephemeris_, "Earth"));
  auto const earth_degrees_of_freedom =
      solar_system_2000_.degrees_of_freedom("Earth");

  Time const sidereal_day = Day * 365.2425 / 366.2425;

  // These data are from https://en.wikipedia.org/wiki/Molniya_orbit.  The
  // eccentricity is from the "External links" section.
  KeplerianElements<ICRS> initial_elements;
  initial_elements.eccentricity = 0.74105;
  initial_elements.mean_motion = 2 * π * Radian / (sidereal_day / 2);
  initial_elements.inclination = ArcSin(2.0 / Sqrt(5.0));
  initial_elements.argument_of_periapsis = -π / 2.0 * Radian;
  initial_elements.longitude_of_ascending_node = 1 * Radian + 60 * Degree;
  initial_elements.mean_anomaly = 2 * Radian;

  MasslessBody const satellite{};
  KeplerOrbit<ICRS> initial_orbit(
      *earth_body, satellite, initial_elements, J2000);
  auto const satellite_state_vectors = initial_orbit.StateVectors(J2000);

  OrbitAnalyser<ICRS> analyser(
      ephemeris_.get(),
      earth_body,
      J2000,
      earth_degrees_of_freedom + satellite_state_vectors);
  analyser.Analyse();
}

TEST_F(OrbitAnalyzerTest, 福爾摩沙衛星二號) {
  auto const earth_body = dynamic_cast_not_null<OblateBody<ICRS> const*>(
      solar_system_2000_.massive_body(*ephemeris_, "Earth"));
  // 1 28254U 04018A   16183.13910299  .00000106  00000-0  10000-3 0  9998
  // 2 28254  98.9332 246.1213 0002410  82.2354 330.6801 14.00730683619455

  // The above two-line elements converted to GCRS state vectors at their epoch
  // with Skyfield.
  // The satellite was operational at that epoch, so we can expect the ground
  // track recurrence to be correctly set up.
  // While the two-line elements can be used to directly obtain accurate periods
  // (Ἰξίων computes κ = 13.999159 from these elements, indicating very good
  // recurrence), their conversion to instantataneous state vectors has
  // uncertainties on the order of a km; indeed, we end up measuring κ = 13.994
  // using the output of Skyfield.
  // As we have no more than 4 sig. dec. on these values, we apply a fudge
  // factor that gets us to a good recurrence, which is the what we want to
  // study with this test.
  // The recurrence is not *too good* either: the fudge factor is a multiple of
  // 1e-5, which corresponds to a precision of 74 mm/s on the velocity; this is
  // approximately the amount of Δv needed to correct for drag over a year on
  // this orbit.
  auto const epoch = "2016-07-01T03:20:18"_UTC;
  RelativeDegreesOfFreedom<ICRS> const satellite_state_vectors = {
      Displacement<ICRS>({-1.73980671e-05 * AstronomicalUnit,
                          -2.43219608e-05 * AstronomicalUnit,
                          3.82473828e-05 * AstronomicalUnit}),
      Velocity<ICRS>({0.00103147 * AstronomicalUnit / Day,
                      0.00327939 * AstronomicalUnit / Day,
                      0.00254849 * AstronomicalUnit / Day}) *
          (1 - 14e-5)};

  ephemeris_->Prolong(epoch);
  auto const earth_degrees_of_freedom =
      ephemeris_->trajectory(earth_body)->EvaluateDegreesOfFreedom(epoch);

  OrbitAnalyser<ICRS> analyser(
      ephemeris_.get(),
      earth_body,
      epoch,
      earth_degrees_of_freedom + satellite_state_vectors);
  analyser.Analyse();
}

TEST_F(OrbitAnalyzerTest, DISABLED_あけぼの) {
  auto const earth_body = dynamic_cast_not_null<OblateBody<ICRS> const*>(
      solar_system_2000_.massive_body(*ephemeris_, "Earth"));
  auto const earth_degrees_of_freedom =
      solar_system_2000_.degrees_of_freedom("Earth");

  // These data are from Michel Capderou (2011), Satellites : de Kepler au GPS,
  // Fig. 9.6 (p.332).
  KeplerianElements<ICRS> initial_elements;
  initial_elements.eccentricity = 0.269151;
  initial_elements.semimajor_axis =  9088.656 * Kilo(Metre);
  initial_elements.inclination = 75.07 * Degree;
  initial_elements.argument_of_periapsis = 149.24 * Degree;
    // TODO(egg): the LAN is probably wrong.
  initial_elements.longitude_of_ascending_node = 98.26 * Degree;
  initial_elements.mean_anomaly = 2 * Radian;

  MasslessBody const satellite{};
  KeplerOrbit<ICRS> initial_orbit(
      *earth_body, satellite, initial_elements, J2000);
  auto const satellite_state_vectors = initial_orbit.StateVectors(J2000);

  OrbitAnalyser<ICRS> analyser(
      ephemeris_.get(),
      earth_body,
      J2000,
      earth_degrees_of_freedom + satellite_state_vectors);
  analyser.Analyse();
}

#endif

}  // namespace astronomy
}  // namespace principia
