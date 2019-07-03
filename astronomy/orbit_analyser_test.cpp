#include "astronomy/orbit_analyser.hpp"

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
using base::make_not_null_unique;
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
using ::testing::Eq;

namespace astronomy {

struct SP3Files {
  std::vector<std::string> names;
  StandardProduct3::Dialect dialect;
};

struct SP3Orbit {
  StandardProduct3::SatelliteIdentifier satellite;
  SP3Files files;
};

std::ostream& operator<<(std::ostream& out, SP3Orbit const& orbit) {
  out << orbit.satellite << " in (";
  for (int i = 0; i < orbit.files.names.size(); ++i) {
    out << orbit.files.names[i] << ", ";
  }
  return out << "interpreted as " << orbit.files.dialect << ")";
}

class OrbitAnalyserTest : public ::testing::TestWithParam<SP3Orbit> {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
    if (ephemeris_ == nullptr) {
      ephemeris_ = solar_system_.MakeEphemeris(
          /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                   /*geopotential_tolerance=*/0x1p-24},
          Ephemeris<ICRS>::FixedStepParameters(
              SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                 Position<ICRS>>(),
              /*step=*/10 * Minute));
      for (char const* const date :
           {"2014-06-01T12:00:00", "2019-01-01T12:00:00"}) {
        LOG(INFO) << "Prolonging ephemeris to " << date << " (TT)";
        ephemeris_->Prolong(ParseTT(date));
      }
    }
  }

  not_null<std::unique_ptr<DiscreteTrajectory<ICRS>>> Trajectory() {
    auto result = make_not_null_unique<DiscreteTrajectory<ICRS>>();
    for (auto const& file : GetParam().files.names) {
      StandardProduct3 sp3(
          SOLUTION_DIR / "astronomy" / "standard_product_3" / file,
          GetParam().files.dialect);
      auto const& orbit = sp3.orbit(GetParam().satellite);
      CHECK_EQ(orbit.size(), 1);
      auto const& arc = *orbit.front();
      for (auto it = arc.Begin(); it != arc.End(); ++it) {
        ephemeris_->Prolong(it.time());
        result->Append(
            it.time(),
            itrs_.FromThisFrameAtTime(it.time())(it.degrees_of_freedom()));
      }
    }
    return result;
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
      sun_,
      J2000,
      earth_trajectory_.EvaluateDegreesOfFreedom(J2000) +
          satellite_state_vectors,
      u8"Молния");
  analyser.Analyse();
}

std::vector<SP3Orbit> const& GNSSOrbits() {
  static const SP3Files files = {{"WUM0MGXFIN_20190970000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20190980000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20190990000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191000000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191010000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191020000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191030000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191040000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191050000_01D_15M_ORB.SP3",
                                  "WUM0MGXFIN_20191060000_01D_15M_ORB.SP3"},
                                 StandardProduct3::Dialect::ChineseMGEX};
  // J07 is not present in the last two files.
  static const SP3Files first_8_files = []() {
    SP3Files result = files;
    result.names.pop_back();
    result.names.pop_back();
    return result;
  }();
  static const std::vector<SP3Orbit> orbits{{
      {{StandardProduct3::SatelliteGroup::北斗, 1}, files},     // GEO 140° E.
      {{StandardProduct3::SatelliteGroup::北斗, 6}, files},     // IGSO 118° E.
      {{StandardProduct3::SatelliteGroup::北斗, 11}, files},    // MEO A07.
      {{StandardProduct3::SatelliteGroup::Galileo, 1}, files},  // MEO B05.
      {{StandardProduct3::SatelliteGroup::GPS, 1}, files},      // MEO D2.
      {{StandardProduct3::SatelliteGroup::みちびき, 1}, files},  // QZO 139° E.
      {{StandardProduct3::SatelliteGroup::みちびき, 7},          // GEO 127° E.
       first_8_files},
      {{StandardProduct3::SatelliteGroup::ГЛОНАСС, 1}, files},  // MEO Plane I.
  }};
  return orbits;
}

std::vector<SP3Orbit> const& SPOT5Orbit() {
  static const SP3Files files = {{"ssasp501.b10170.e10181.D__.sp3"},
                                 StandardProduct3::Dialect::Standard};
  static const std::vector<SP3Orbit> orbits{{
      {{StandardProduct3::SatelliteGroup::General, 94}, files},
  }};
  return orbits;
}

std::vector<SP3Orbit> const& Sentinel3AOrbit() {
  static const SP3Files files = {{"ssas3a20.b18358.e19003.DG_.sp3"},
                                 StandardProduct3::Dialect::Standard};
  static const std::vector<SP3Orbit> orbits{{
      {{StandardProduct3::SatelliteGroup::General, 74}, files},
  }};
  return orbits;
}

std::vector<SP3Orbit> const& TOPEXPoséidonOrbit() {
  static const SP3Files files = {{"grgtop03.b97344.e97348.D_S.sp3"},
                                 StandardProduct3::Dialect::GRGS};
  static const std::vector<SP3Orbit> orbits{{
      {{StandardProduct3::SatelliteGroup::General, 1}, files},
  }};
  return orbits;
}

INSTANTIATE_TEST_CASE_P(GNSS,
                        OrbitAnalyserTest,
                        ::testing::ValuesIn(GNSSOrbits()));

INSTANTIATE_TEST_CASE_P(SPOT5,
                        OrbitAnalyserTest,
                        ::testing::ValuesIn(SPOT5Orbit()));

INSTANTIATE_TEST_CASE_P(Sentinel3A,
                        OrbitAnalyserTest,
                        ::testing::ValuesIn(Sentinel3AOrbit()));

INSTANTIATE_TEST_CASE_P(TOPEXPoséidon,
                        OrbitAnalyserTest,
                        ::testing::ValuesIn(TOPEXPoséidonOrbit()));

TEST_P(OrbitAnalyserTest, DoTheAnalysis) {
  OrbitAnalyser<ICRS> analyser(
      ephemeris_.get(),
      earth_,
      sun_,
      *Trajectory(),
      (std::stringstream() << GetParam().satellite).str());
  analyser.RecomputeProperties();
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
      sun_,
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
      earth_, sun_,
      epoch,
      earth_trajectory_.EvaluateDegreesOfFreedom(epoch) +
          satellite_state_vectors,
      "FormoSat-2");
  analyser.Analyse();
}

#endif

}  // namespace astronomy
}  // namespace principia
