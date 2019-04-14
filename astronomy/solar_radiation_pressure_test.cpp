
#include <limits>

#include "absl/strings/str_cat.h"
#include "astronomy/standard_product_3.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"

namespace principia {
namespace astronomy {

using base::not_null;
using geometry::Instant;
using geometry::Position;
using geometry::Vector;
using integrators::EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Fine1987RKNG34;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodySurfaceDynamicFrame;
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::Acceleration;
using quantities::Length;
using quantities::Square;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;

using ::testing::TestWithParam;

namespace {

struct StandardProduct3Args {
  std::filesystem::path filename;
  StandardProduct3::Dialect dialect = StandardProduct3::Dialect::Standard;
};

std::ostream& operator<<(std::ostream& out, StandardProduct3Args const& args) {
  return out << args.filename << " interpreted as " << args.dialect;
}

}  // namespace

class SolarRadiationPressureTest : public TestWithParam<StandardProduct3Args> {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
    if (static_ephemeris_ == nullptr) {
      static_ephemeris_ = solar_system_2010_.MakeEphemeris(
          /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                   /*geopotential_tolerance=*/0x1p-24},
          Ephemeris<ICRS>::FixedStepParameters(
              SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                 Position<ICRS>>(),
              /*step=*/10 * Minute));
      for (int year = 2011; year <= 2020; ++year) {
        std::string date = absl::StrCat(year, "-01-01T00:00:00");
        LOG(INFO) << "Prolonging to " << date << " (UTC)";
        static_ephemeris_->Prolong(ParseUTC(date));
      }
    }
  }

 protected:
  SolarRadiationPressureTest()
      : sp3_(GetParam().filename, GetParam().dialect),
        ephemeris_(*static_ephemeris_),
        earth_(solar_system_2010_.rotating_body(ephemeris_, "Earth")),
        sun_(solar_system_2010_.rotating_body(ephemeris_, "Sun")),
        earth_trajectory_(*ephemeris_.trajectory(earth_)),
        sun_trajectory_(*ephemeris_.trajectory(sun_)),
        itrs_(&ephemeris_, earth_) {}

  Square<Length> Residual(
      StandardProduct3::SatelliteIdentifier const& satellite,
      Ephemeris<ICRS>::GeneralizedIntrinsicAcceleration const&
          solar_radiation_pressure) const {
    DiscreteTrajectory<ICRS> integrated;
    integrated.Append(
        sp3_.orbit(satellite).front()->Begin().time(),
        itrs_.FromThisFrameAtTime(
            sp3_.orbit(satellite).front()->Begin().time())(
            sp3_.orbit(satellite).front()->Begin().degrees_of_freedom()));
    Square<Length> residual;
    for (not_null<DiscreteTrajectory<ITRS> const*> const arc :
         sp3_.orbit(satellite)) {
      for (auto it = arc->Begin(); it != arc->End(); ++it) {
        Ephemeris<ICRS>::GeneralizedAdaptiveStepParameters parameters(
            EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                Fine1987RKNG34,
                Position<ICRS>>(),
            /*max_steps=*/std::numeric_limits<int64_t>::max(),
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Micro(Metre) / Second);
        ephemeris_.FlowWithAdaptiveStep(
            &integrated,
            solar_radiation_pressure,
            it.time(),
            parameters,
            /*max_ephemeris_steps=*/std::numeric_limits<int64_t>::max(),
            /*last_point_only=*/true);
        residual += (itrs_.ToThisFrameAtTime(it.time())(
                         integrated.last().degrees_of_freedom()).position() -
                     it.degrees_of_freedom().position()).Norm²();
      }
    }
    return residual;
  }

 private:
  static SolarSystem<ICRS> const solar_system_2010_;
  static std::unique_ptr<Ephemeris<ICRS>> static_ephemeris_;

 protected:
  StandardProduct3 const sp3_;
  Ephemeris<ICRS>& ephemeris_;
  not_null<RotatingBody<ICRS> const*> const earth_;
  not_null<RotatingBody<ICRS> const*> const sun_;
  ContinuousTrajectory<ICRS> const& earth_trajectory_;
  ContinuousTrajectory<ICRS> const& sun_trajectory_;
  BodySurfaceDynamicFrame<ICRS, ITRS> const itrs_;
};

SolarSystem<ICRS> const SolarRadiationPressureTest::solar_system_2010_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2455200_500000000.proto.txt");
std::unique_ptr<Ephemeris<ICRS>> SolarRadiationPressureTest::static_ephemeris_;

// Optimize a constant radial force.
TEST_P(SolarRadiationPressureTest, RadialForce) {
  for (auto const& satellite : sp3_.satellites()) {
    Acceleration radial_acceleration = 1e-7 * Metre / Pow<2>(Second);
    auto solar_radiation_pressure =
        [](Instant const& t,
           DegreesOfFreedom<ICRS> const& dof) -> Vector<Acceleration, ICRS> {
    };
  }
}

}  // namespace astronomy
}  // namespace principia
