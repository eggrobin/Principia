#include "astronomy/orbit_analyser.hpp"

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
using base::not_null;
using base::OFStream;
using geometry::Displacement;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodySurfaceDynamicFrame;
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
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
using quantities::si::Deci;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;

namespace astronomy {

class OrbitAnalyserTest : public ::testing::Test {
 protected:
  OrbitAnalyserTest()
      : earth_(dynamic_cast_not_null<OblateBody<ICRS> const*>(
            solar_system_.massive_body(*ephemeris_, "Earth"))),
        earth_trajectory_(*ephemeris_->trajectory(earth_)),
        itrs_(ephemeris_.get(), earth_) {}

  static void SetUpTestCase() {
    google::LogToStderr();
    ephemeris_ = solar_system_.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<ICRS>>(),
            /*step=*/10 * Minute));
  }

  not_null<OblateBody<ICRS> const*> const earth_;
  ContinuousTrajectory<ICRS> const& earth_trajectory_;
  BodySurfaceDynamicFrame<ICRS, ITRS> itrs_;

  static SolarSystem<ICRS> solar_system_;
  static std::unique_ptr<Ephemeris<ICRS>> ephemeris_;
};

SolarSystem<ICRS> OrbitAnalyserTest::solar_system_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2436145_604166667.proto.txt");
std::unique_ptr<Ephemeris<ICRS>> OrbitAnalyserTest::ephemeris_;

#if !defined(_DEBUG)

TEST_F(OrbitAnalyserTest, DISABLED_Молния) {
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
      *earth_, satellite, initial_elements, J2000);
  auto const satellite_state_vectors = initial_orbit.StateVectors(J2000);

  OrbitAnalyser<ICRS> analyser(
      ephemeris_.get(),
      earth_,
      J2000,
      earth_trajectory_.EvaluateDegreesOfFreedom(J2000) +
          satellite_state_vectors);
  analyser.Analyse();
}

TEST_F(OrbitAnalyserTest, TOPEXPoséidon) {
  // Initial state from DORIS products.
  // https://ids-doris.org/ids/organization/data-centers.html.
  // See also the definition of the SP3 format
  // ftp://igs.org/pub/data/format/sp3c.txt.

  // grgtop03.b97344.e97348.D_S.sp3, header and first record, from
  // ftp://doris.ign.fr/pub/doris/products/orbits/grg/top/.
  // Note: LCA stands for LEGOS/CLS AC, where LEGOS is Laboratoire d’Études en
  // Géophysique et Océanographie Spatiales, CLS is Collecte Localisation
  // Satellites, and AC is Analysis Centre.
  // However, CNES/CLS transitioned their AC acronym from LCA to GRG in 2014
  // (where GRG probably stands for Groupe de Recherche de Géodesie Spatiale),
  // hence the mismatch between the file name and the agency field of the SP3
  // header below.
  // #cV1997 12 10 12  0  0.00000000    5046 DORIS ITR05 FIT  LCA
  // ##  935 302400.00000000    60.00000000 50792 0.5000000000000
  // +    1   L01  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  // [... Lines 4 through 12 omitted                         ...]
  // %c L  cc TAI ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  // %c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  // %f  0.0000000  0.000000000  0.00000000000  0.000000000000000
  // [... Lines 15 through 18 omitted                        ...]
  // /* CNES/CLS - TOULOUSE, FRANCE ccccc ccccc ccccc ccccc ccccc
  // /* CCCCCCCCCCCCCC Contact : Laurent SOUDARIN CCCCCCCCCCCCCCC
  // /* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  // /* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  // *  1997 12 10 12  0  0.00000000
  // PL01  -3091.510103   1090.750605  -6985.258847 999999.999999
  // VL01   -414.743700  -6885.134975   -890.998602 999999.999999

  constexpr Instant initial_time = "1997-12-10T12:00:00,000"_TAI;
  DegreesOfFreedom<ITRS> const initial_dof = {
      ITRS::origin + Displacement<ITRS>({-3091.510103 * Kilo(Metre),
                                          1090.750605 * Kilo(Metre),
                                         -6985.258847 * Kilo(Metre)}),
      Velocity<ITRS>({ -414.743700 * Deci(Metre) / Second,
                      -6885.134975 * Deci(Metre) / Second,
                       -890.998602 * Deci(Metre) / Second})};

  ephemeris_->Prolong(initial_time);
  OrbitAnalyser<ICRS> analyser(
      ephemeris_.get(),
      earth_,
      initial_time,
      itrs_.FromThisFrameAtTime(initial_time)(initial_dof));
  analyser.Analyse();
}

TEST_F(OrbitAnalyserTest, 福爾摩沙衛星二號) {
  // 1 28254U 04018A   16183.13910299  .00000106  00000-0  10000-3 0  9998
  // 2 28254  98.9332 246.1213 0002410  82.2354 330.6801 14.00730683619455

  // The above two-line elements converted to GCRS state vectors at their epoch
  // with Skyfield.
  // The satellite was operational at that epoch, so we can expect the ground
  // track recurrence to be correctly set up.
  // While the two-line elements can be used to directly obtain accurate periods
  // (Ἰξίων computes κ = 13.999159 from these elements, indicating very good
  // recurrence), their conversion to instantataneous state vectors has
  // uncertainties on the order of a kilometre; indeed, we end up measuring
  // κ = 13.994 using the output of Skyfield.
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

  OrbitAnalyser<ICRS> analyser(
      ephemeris_.get(),
      earth_,
      epoch,
      earth_trajectory_.EvaluateDegreesOfFreedom(epoch) +
          satellite_state_vectors);
  analyser.Analyse();
}

#endif

}  // namespace astronomy
}  // namespace principia
