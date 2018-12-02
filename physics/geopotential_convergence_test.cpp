﻿
#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/statistics.hpp"

namespace principia {

using astronomy::ICRS;
using astronomy::J2000;
using base::dynamic_cast_not_null;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order12;
using quantities::ArcSin;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::EarthEquatorialRadius;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;

namespace physics {

class GeopotentialConvergenceTest : public ::testing::Test {
 protected:
  static void SetUpTestCase() {
    google::LogToStderr();
    ephemeris_ = solar_system_2000_.MakeEphemeris(
        Ephemeris<ICRS>::AccuracyParameters(
            /*fitting_tolerance=*/5 * Milli(Metre),
            /*geopotential_tolerance=*/0),
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<ICRS>>(),
            /*step=*/10 * Minute));
  }

  static SolarSystem<ICRS> solar_system_2000_;
  static std::unique_ptr<Ephemeris<ICRS>> ephemeris_;
};

SolarSystem<ICRS> GeopotentialConvergenceTest::solar_system_2000_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2451545_000000000.proto.txt");
std::unique_ptr<Ephemeris<ICRS>> GeopotentialConvergenceTest::ephemeris_;


TEST_F(GeopotentialConvergenceTest, EccentricOrbit) {
  auto const earth_body = dynamic_cast_not_null<OblateBody<ICRS> const*>(
      solar_system_2000_.massive_body(*ephemeris_, "Earth"));
  auto const earth_degrees_of_freedom =
      solar_system_2000_.degrees_of_freedom("Earth");
  auto const& geopotential = ephemeris_->geopotential(earth_body);

  Time const sidereal_day = Day * 365.2425 / 366.2425;

  // An orbit with the Молния inclination, but with a much longer period, so as
  // to cross the outer sectoral threshold.
  KeplerianElements<ICRS> initial_elements;
  initial_elements.periapsis_distance = EarthEquatorialRadius + 1000 * Kilo(Metre);
  initial_elements.mean_motion = 2.0 * π * Radian / (10 * sidereal_day);
  initial_elements.inclination = ArcSin(2.0 / Sqrt(5.0));
  initial_elements.argument_of_periapsis = -π / 2.0 * Radian;
  initial_elements.longitude_of_ascending_node = 1 * Radian;
  initial_elements.mean_anomaly = 2 * Radian;

  MasslessBody const satellite{};
  KeplerOrbit<ICRS> initial_orbit(
      *earth_body, satellite, initial_elements, J2000);
  auto const satellite_state_vectors = initial_orbit.StateVectors(J2000);

  EXPECT_THAT(*initial_orbit.elements_at_epoch().apoapsis_distance,
              AllOf(Gt(geopotential.sectoral_damping().outer_threshold()),
                    Lt(geopotential.degree_damping()[2].inner_threshold())));
  EXPECT_THAT(*initial_orbit.elements_at_epoch().periapsis_distance,
              Lt(geopotential.degree_damping().back().inner_threshold()));

  auto flow_orbit_1_year = [this,
                            &earth_degrees_of_freedom,
                            &satellite_state_vectors](
                               Time const integration_step) {
    Time const integration_duration = 1.0 * JulianYear;
    DiscreteTrajectory<ICRS> trajectory;
    trajectory.Append(J2000,
                      earth_degrees_of_freedom + satellite_state_vectors);
    auto const instance = ephemeris_->NewInstance(
        {&trajectory},
        Ephemeris<ICRS>::NoIntrinsicAccelerations,
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                               Position<ICRS>>(),
            integration_step));

    ephemeris_->FlowWithFixedStep(J2000 + integration_duration, *instance);
    return trajectory.last().degrees_of_freedom().position();
  };

  Position<ICRS> const reference_solution = flow_orbit_1_year(1 * Second);
  double const increment = std::pow(2, 0x1p-8);
  for (Time Δt = 50 * Second; Δt < 100 * Second; Δt *= increment) {
    LOG(ERROR) << Δt << ": "
                     << (flow_orbit_1_year(Δt) - reference_solution).Norm();
  }
}

}  // namespace physics
}  // namespace principia
