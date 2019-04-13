﻿#include "astronomy/orbit_analyser.hpp"

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/standard_product_3.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
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
using quantities::si::Micro;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using ::testing::Eq;

namespace astronomy {

struct SP3Orbit {
  std::string_view filename;
  StandardProduct3::SatelliteIdentifier satellite;
  StandardProduct3::Dialect dialect = StandardProduct3::Dialect::Standard;
};

class OrbitAnalyserTest : public ::testing::TestWithParam<SP3Orbit> {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
    if (ephemeris_ == nullptr) {
      ephemeris_ = solar_system_.MakeEphemeris(
          /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Micro(Metre),
                                   /*geopotential_tolerance=*/0x1p-24},
          Ephemeris<ICRS>::FixedStepParameters(
              SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                 Position<ICRS>>(),
              /*step=*/1 * Minute));
      for (char const* const date :
           {"2014-06-01T12:00:00", "2019-01-01T12:00:00"}) {
        LOG(INFO) << "Prolonging ephemeris to " << date << " (TT)";
        ephemeris_->Prolong(ParseTT(date));
      }
    }
  }

 protected:
  OrbitAnalyserTest()
      : earth_(dynamic_cast_not_null<OblateBody<ICRS> const*>(
            solar_system_.massive_body(*ephemeris_, "Earth"))),
        sun_(dynamic_cast_not_null<OblateBody<ICRS> const*>(
            solar_system_.massive_body(*ephemeris_, "Sun"))),
        earth_trajectory_(*ephemeris_->trajectory(earth_)),
        itrs_(ephemeris_.get(), earth_) {}

  not_null<OblateBody<ICRS> const*> const earth_;
  not_null<OblateBody<ICRS> const*> const sun_;
  ContinuousTrajectory<ICRS> const& earth_trajectory_;
  BodySurfaceDynamicFrame<ICRS, ITRS> itrs_;

  static SolarSystem<ICRS> solar_system_;
  static std::unique_ptr<Ephemeris<ICRS>> ephemeris_;
};

SolarSystem<ICRS> OrbitAnalyserTest::solar_system_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2455200_500000000.proto.txt");
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
          satellite_state_vectors,
      u8"Молния");
  analyser.Analyse();
}

INSTANTIATE_TEST_CASE_P(LAGEOS2,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "ilrsa.orb.lageos2.160319.v35.sp3",
                            {StandardProduct3::SatelliteGroup::General, 52},
                            StandardProduct3::Dialect::ILRSA}));
INSTANTIATE_TEST_CASE_P(GPS,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "nga20342.eph",
                            {StandardProduct3::SatelliteGroup::GPS, 1}}));
INSTANTIATE_TEST_CASE_P(Galileo,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "COD0MGXFIN_20183640000_01D_05M_ORB.SP3",
                            {StandardProduct3::SatelliteGroup::Galileo, 1}}));
INSTANTIATE_TEST_CASE_P(ГЛОНАСС,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "COD0MGXFIN_20183640000_01D_05M_ORB.SP3",
                            {StandardProduct3::SatelliteGroup::ГЛОНАСС, 1}}));
// Whereas the AIUB, GFZ, GRGS, SHAO, TUM, and WHU MGEX analysis centres all
// provide orbits for 北斗 satellites, only GFZ, SHAO, TUM, and WHU provide orbits
// for the geostationary satellites C01 ‥ C05.  Of those, only WHU provides
// final analysis products; GFZ, SHAO, and TUM provide rapid analysis products.
INSTANTIATE_TEST_CASE_P(北斗Geostationary80,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "WUM0MGXFIN_20190270000_01D_15M_ORB.SP3",
                            {StandardProduct3::SatelliteGroup::北斗, 2},
                            StandardProduct3::Dialect::ChineseMGEX}));
INSTANTIATE_TEST_CASE_P(北斗Geostationary160,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "WUM0MGXFIN_20190270000_01D_15M_ORB.SP3",
                            {StandardProduct3::SatelliteGroup::北斗, 4},
                            StandardProduct3::Dialect::ChineseMGEX}));
INSTANTIATE_TEST_CASE_P(北斗InclinedGeosynchronous,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "COD0MGXFIN_20183640000_01D_05M_ORB.SP3",
                            {StandardProduct3::SatelliteGroup::北斗, 6}}));
INSTANTIATE_TEST_CASE_P(北斗MediumEarthOrbit,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "COD0MGXFIN_20183640000_01D_05M_ORB.SP3",
                            {StandardProduct3::SatelliteGroup::北斗, 11}}));
INSTANTIATE_TEST_CASE_P(QZSS,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "COD0MGXFIN_20183640000_01D_05M_ORB.SP3",
                            {StandardProduct3::SatelliteGroup::みちびき, 1}}));
// Whereas JAXA and the AIUB, GFZ, TUM, and WHU MGEX analysis centres all
// provide orbits for QZSS satellites, only GFZ and WHU provide orbits for the
// geostationary satellite J07.  Of those, only WHU provides final analysis
// products; GFZ provides rapid analysis products.
INSTANTIATE_TEST_CASE_P(QZSSGeostationary,
                        OrbitAnalyserTest,
                        ::testing::Values(SP3Orbit{
                                "WUM0MGXFIN_20190270000_01D_15M_ORB.SP3",
                            {StandardProduct3::SatelliteGroup::みちびき, 7},
                            StandardProduct3::Dialect::ChineseMGEX}));

TEST_P(OrbitAnalyserTest, Residuals) {
  StandardProduct3 sp3(
      SOLUTION_DIR / "astronomy" / "standard_product_3" / GetParam().filename,
      GetParam().dialect);
  StandardProduct3::SatelliteIdentifier const& satellite = GetParam().satellite;
  std::string name = (std::stringstream() << satellite).str();

  std::vector<Displacement<ICRS>> icrs_residuals;
  std::vector<Displacement<ITRS>> itrs_residuals;
  std::vector<Length> altitudes;
  std::vector<Angle> sun_earth_satellite_angle;
  std::vector<double> times;
  DiscreteTrajectory<ICRS> trajectory;
  auto const start_of_first_arc = sp3.orbit(satellite).front()->Begin();
  ephemeris_->Prolong(start_of_first_arc.time());
  trajectory.Append(start_of_first_arc.time(),
                    itrs_.FromThisFrameAtTime(start_of_first_arc.time())(
                        start_of_first_arc.degrees_of_freedom()));

  for (auto const& arc : sp3.orbit(satellite)) {
    for (auto it = arc->Begin();;) {
      ephemeris_->Prolong(it.time());

      if (++it == arc->End()) {
        break;
      }

      Ephemeris<ICRS>::AdaptiveStepParameters parameters(
          integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator<
             integrators::methods::DormandالمكاوىPrince1986RKN434FM,
              Position<ICRS>>(),
          /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
          1 * Micro(Metre),
          1 * Micro(Metre) / Second);
      ephemeris_->FlowWithAdaptiveStep(
          &trajectory,
          Ephemeris<ICRS>::NoIntrinsicAcceleration,
          it.time(),
          parameters,
          /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
          /*last_point_only=*/true);

      times.push_back((it.time() - J2000) / Day);
      icrs_residuals.push_back(
          trajectory.last().degrees_of_freedom().position() -
          itrs_.FromThisFrameAtTime(it.time())(it.degrees_of_freedom())
              .position());
      itrs_residuals.push_back(itrs_
                                   .ToThisFrameAtTime(it.time())(
                                       trajectory.last().degrees_of_freedom())
                                   .position() -
                               it.degrees_of_freedom().position());
      sun_earth_satellite_angle.push_back(geometry::AngleBetween(
          itrs_.FromThisFrameAtTime(it.time())(it.degrees_of_freedom())
                  .position() -
              ephemeris_->trajectory(earth_)->EvaluatePosition(it.time()),
          ephemeris_->trajectory(sun_)->EvaluatePosition(it.time()) -
              ephemeris_->trajectory(earth_)->EvaluatePosition(it.time())));
      altitudes.push_back(
          (itrs_.FromThisFrameAtTime(it.time())(it.degrees_of_freedom())
               .position() -
           ephemeris_->trajectory(earth_)->EvaluatePosition(it.time()))
              .Norm() -
          earth_->mean_radius());
    }
    base::OFStream f(SOLUTION_DIR / ("residuals_" + name));
    f << mathematica::Assign(
        mathematica::Apply("residualsICRS", {mathematica::Escape(name)}),
        mathematica::Apply(
            "Transpose",
            {mathematica::Apply(
                "List",
                {mathematica::ToMathematica(times),
                 mathematica::ToMathematica(
                     mathematica::ExpressIn(Metre, icrs_residuals))})}));
    f << mathematica::Assign(
        mathematica::Apply("residualsITRS", {mathematica::Escape(name)}),
        mathematica::Apply(
            "Transpose",
            {mathematica::Apply(
                "List",
                {mathematica::ToMathematica(times),
                 mathematica::ToMathematica(
                     mathematica::ExpressIn(Metre, itrs_residuals))})}));
    f << mathematica::Assign(
        mathematica::Apply("sunEarthSatAngles", {mathematica::Escape(name)}),
        mathematica::Apply(
            "Transpose",
            {mathematica::Apply(
                "List",
                {mathematica::ToMathematica(times),
                 mathematica::ToMathematica(
                     mathematica::ExpressIn(Radian, sun_earth_satellite_angle))})}));
    f << mathematica::Assign(
        mathematica::Apply("altitudes", {mathematica::Escape(name)}),
        mathematica::Apply(
            "Transpose",
            {mathematica::Apply(
                "List",
                {mathematica::ToMathematica(times),
                 mathematica::ToMathematica(
                     mathematica::ExpressIn(Kilo(Metre), altitudes))})}));
  }
}

TEST_P(OrbitAnalyserTest, GNSS) {
  StandardProduct3 sp3(
      SOLUTION_DIR / "astronomy" / "standard_product_3" / GetParam().filename,
      GetParam().dialect);
  StandardProduct3::SatelliteIdentifier const& satellite = GetParam().satellite;

  auto it = sp3.orbit(satellite).front()->Begin();
  for (int i = 0; i < 10; ++i) {
    ++it;
  }

  ephemeris_->Prolong(it.time());
  OrbitAnalyser<ICRS> analyser(
      ephemeris_.get(),
      earth_,
      it.time(),
      itrs_.FromThisFrameAtTime(it.time())(it.degrees_of_freedom()),
      (std::stringstream() << satellite).str());
  analyser.Analyse();
}

TEST_F(OrbitAnalyserTest, TOPEXPoséidon) {
  StandardProduct3 sp3(SOLUTION_DIR / "astronomy" / "standard_product_3" /
                           "grgtop03.b97344.e97348.D_S.sp3",
                       StandardProduct3::Dialect::GRGS);
  StandardProduct3::SatelliteIdentifier topex_poséidon{
      StandardProduct3::SatelliteGroup::General, 1};

  Instant const initial_time =
      sp3.orbit(topex_poséidon).front()->Begin().time();

  EXPECT_THAT(initial_time, Eq("1997-12-10T12:00:00,000"_TAI));


  ephemeris_->Prolong(initial_time);

  OrbitAnalyser<ICRS> analyser(
      ephemeris_.get(),
      earth_,
      initial_time,
      itrs_.FromThisFrameAtTime(initial_time)(
          sp3.orbit(topex_poséidon).front()->Begin().degrees_of_freedom()),
      "TOPEX/Poseidon");
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
          satellite_state_vectors,
      "FormoSat-2");
  analyser.Analyse();
}

#endif

}  // namespace astronomy
}  // namespace principia
