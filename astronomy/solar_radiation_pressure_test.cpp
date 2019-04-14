
#include <limits>
#include <random>

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
using geometry::AngleBetween;
using geometry::Bivector;
using geometry::Displacement;
using geometry::Instant;
using geometry::Normalize;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using geometry::Wedge;
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
using quantities::ArcSin;
using quantities::Length;
using quantities::Pow;
using quantities::Sqrt;
using quantities::Square;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;

using ::testing::TestWithParam;
using ::testing::ValuesIn;

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
      for (int year = 2010; year <= 2019; ++year) {
        std::string date = absl::StrCat(year, "-06-01T00:00:00");
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

  Length Residual(
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
    int n = 0;
    for (not_null<DiscreteTrajectory<ITRS> const*> const arc :
         sp3_.orbit(satellite)) {
      for (auto it = arc->Begin(); it != arc->End(); ++it, ++n) {
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
        residual = std::max(
            residual,
            (itrs_.ToThisFrameAtTime(it.time())(
                  integrated.last().degrees_of_freedom()).position() -
             it.degrees_of_freedom().position()).Norm²());
      }
    }
    return Sqrt(residual);
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

INSTANTIATE_TEST_CASE_P(WHUCODE,
                        SolarRadiationPressureTest,
                        ValuesIn(std::vector<StandardProduct3Args>{
                            {SOLUTION_DIR / "astronomy" / "standard_product_3" /
                                 "WUM0MGXFIN_20190270000_01D_15M_ORB.SP3",
                             StandardProduct3::Dialect::ChineseMGEX},
                            {SOLUTION_DIR / "astronomy" / "standard_product_3" /
                             "COD0MGXFIN_20181260000_01D_05M_ORB.SP3"},
                        }));

// Optimize a constant radial force.
TEST_P(SolarRadiationPressureTest, RadialForce) {
  for (auto const& satellite : sp3_.satellites()) {
    struct DYB;
    Vector<Acceleration, DYB> candidate_radial_acceleration;
    std::mt19937_64 engine(1729);
    double temperature = 1;
    auto solar_radiation_pressure =
        [this, &candidate_radial_acceleration](
            Instant const& t, DegreesOfFreedom<ICRS> const& satellite_dof)
        -> Vector<Acceleration, ICRS> {
      Position<ICRS> const earth_position = earth_trajectory_.EvaluatePosition(t);
      Position<ICRS> const sun_position = sun_trajectory_.EvaluatePosition(t);
      Displacement<ICRS> const satellite_to_sun =
          sun_position - satellite_dof.position();
      // Notation from Arnold et al. (2015), CODE’s new solar radiation pressure
      // model for GNSS orbit determination.
      Displacement<ICRS> const r = satellite_dof.position() - earth_position;
      Vector<double, ICRS> const eD = Normalize(satellite_to_sun);
      Vector<double, ICRS> const er = Normalize(r);
      Bivector<double, ICRS> const eY = -Normalize(Wedge(er, eD));
      Vector<double, ICRS> const eB = eD * eY;
      Rotation<DYB, ICRS> to_icrs(eD, eY, eB);
      // Shadow of a point sun; alternatively, we could rigorously compute the
      // umbra and penumbra.
      if (AngleBetween(satellite_to_sun, -r) <
          ArcSin(earth_->mean_radius() / r.Norm())) {
        return Vector<Acceleration, ICRS>();
      }
      return to_icrs(candidate_radial_acceleration);
    };
    LOG(INFO) << satellite << " no radial acceleration: "
              << Residual(satellite, solar_radiation_pressure);
    Vector<Acceleration, DYB> radial_acceleration(
        {-1e-7 * Metre / Pow<2>(Second),
         0 * Metre / Pow<2>(Second),
         0 * Metre / Pow<2>(Second)});
    candidate_radial_acceleration = radial_acceleration;
    Length last_residual = Residual(satellite, solar_radiation_pressure);

    for (;;) {
      candidate_radial_acceleration =
          radial_acceleration +
          Vector<Acceleration, DYB>({std::normal_distribution()(engine) *
                                         (5e-9 * Metre / Pow<2>(Second)),
                                     std::normal_distribution()(engine) *
                                         (1e-9 * Metre / Pow<2>(Second)),
                                     std::normal_distribution()(engine) *
                                         (1e-9 * Metre / Pow<2>(Second))});
      Length const candidate_residual =
          Residual(satellite, solar_radiation_pressure);
      if (candidate_residual < last_residual ||
          std::uniform_real_distribution()(engine) < temperature) {
        LOG(INFO) << satellite << " " << candidate_residual << " "
                  << candidate_radial_acceleration << " T = " << temperature;
        last_residual = candidate_residual;
        radial_acceleration = candidate_radial_acceleration;
      }
      temperature /= 1 + temperature / 10;
      if (temperature < 0.01) {
        break;
      }
    }
  }
}

}  // namespace astronomy
}  // namespace principia
