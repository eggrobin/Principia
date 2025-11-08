#include "astronomy/frames.hpp"
#include "base/bundle.hpp"
#include "astronomy/orbit_ground_track.hpp"
#include "astronomy/orbital_elements.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "physics/solar_system.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "base/status_utilities.hpp"
#include "mathematica/logger.hpp"
#include "physics/apsides.hpp"
#include "numerics/root_finders.hpp"
#include "numerics/gradient_descent.hpp"
#include "numerics/global_optimization.hpp"

namespace principia {
namespace astronomy {

using namespace principia::quantities::_arithmetic;
using namespace principia::base::_bundle;
using namespace principia::astronomy::_orbital_elements;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_time_scales;
using namespace principia::geometry::_interval;
using namespace principia::geometry::_space;
using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_orbit_ground_track;
using namespace principia::base::_not_null;
using namespace principia::integrators::_methods;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_constants;
using namespace principia::quantities::_si;
using namespace principia::geometry::_frame;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_astronomy;
using namespace principia::numerics::_elementary_functions;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;
using namespace principia::integrators::_embedded_explicit_generalized_runge_kutta_nyström_integrator;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::geometry::_grassmann;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_trajectory;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_body_centred_body_direction_reference_frame;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_solar_system;
using namespace principia::geometry::_instant;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::numerics::_root_finders;
using namespace principia::physics::_apsides;
using namespace principia::numerics::_gradient_descent;

class FishyTest : public ::testing::Test {
 protected:
  FishyTest() {
    google::LogToStderr();
  }

  using LVLH = Frame<struct LVLHTag, Arbitrary, Handedness::Right>;

};

std::string PrettyMeanSolarTimeInterval(Interval<Angle> const& times) {
  Time midpoint_remainder = times.midpoint() * (24 * Hour / (2 * π * Radian));
  Time half_measure_remainder =
      times.measure() / 2 * (24 * Hour / (2 * π * Radian));
  double midpoint_h = std::floor(midpoint_remainder / Hour);
  midpoint_remainder -= midpoint_h * Hour;
  double midpoint_m = std::floor(midpoint_remainder / Minute);
  midpoint_remainder -= midpoint_m * Minute;
  double midpoint_s = midpoint_remainder / Second;
  double half_measure_h = std::floor(half_measure_remainder / Hour);
  half_measure_remainder -= half_measure_h * Hour;
  double half_measure_m = std::floor(half_measure_remainder / Minute);
  half_measure_remainder -= half_measure_m * Minute;
  double half_measure_s = half_measure_remainder / Second;
  return std::format("{:02}:{:02}:{:02.3f}±{:02}:{:02}:{:02.3f}",
                     midpoint_h,
                     midpoint_m,
                     midpoint_s,
                     half_measure_h,
                     half_measure_m,
                     half_measure_s);
}

struct FishyParameters {
  std::string name;
  int max_degree;
  bool zonal_only;
  bool srp;
  Instant min_t0;
  double aspect_ratio;
};

TEST_F(FishyTest, FishyConstant) {
  SolarSystem<ICRS> solar_system_1950_(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  for (auto const& body : solar_system_1950_.MakeAllMassiveBodies()) {
    if (body->name() != "Earth") {
      solar_system_1950_.RemoveMassiveBody(body->name());
    }
  }
  solar_system_1950_.LimitOblatenessToDegree("Earth", 2);
  solar_system_1950_.LimitOblatenessToZonal("Earth");
  Logger logger(TEMP_DIR / "fischbacher_field.wl", /*make_unique=*/false);

  not_null<std::unique_ptr<Ephemeris<ICRS>>> ephemeris_ =
      solar_system_1950_.MakeEphemeris(
          /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                   /*geopotential_tolerance=*/0x1p-24},
          Ephemeris<ICRS>::FixedStepParameters(
              SymmetricLinearMultistepIntegrator<
                  QuinlanTremaine1990Order12,
                  Ephemeris<ICRS>::NewtonianMotionEquation>(),
              /*step=*/10 * Minute));
  RotatingBody<ICRS> const& earth_ =
      *solar_system_1950_.rotating_body(*ephemeris_, "Earth");
  BodyCentredNonRotatingReferenceFrame<ICRS, GCRS> const gcrs_(ephemeris_.get(),
                                                               &earth_);
  Instant const t_max = J2000 + 2 * Day;
  CHECK_OK(ephemeris_->Prolong(t_max));
  double inclination = 97.95;
  double lan = 10;
    std::vector<DiscreteTrajectory<ICRS>> icrs_trajectories;
    icrs_trajectories.reserve(1600);
    icrs_trajectories.emplace_back();
    auto& central_icrs_trajectory = icrs_trajectories.back();
    KeplerOrbit<ICRS> initial_osculating_orbit(
        earth_,
        MasslessBody{},
        KeplerianElements<ICRS>{
            .eccentricity = 0,
            .semimajor_axis = earth_.mean_radius() + 650 * Kilo(Metre),
            .inclination = inclination * Degree,
            .longitude_of_ascending_node = lan * Degree,
            .argument_of_periapsis = 0 * Degree,
            .mean_anomaly = 0 * Degree,
        },
        J2000);
    CHECK_OK(central_icrs_trajectory.Append(
        J2000,
        ephemeris_->trajectory(&earth_)->EvaluateDegreesOfFreedom(J2000) +
            initial_osculating_orbit.StateVectors(J2000)));
    BodyCentredBodyDirectionReferenceFrame<ICRS, LVLH> lvlh(
        ephemeris_.get(),
        [&]() -> Trajectory<ICRS> const& { return central_icrs_trajectory; },
        &earth_);

    Instant const t0 = J2000;

    
    CHECK_OK(ephemeris_->FlowWithAdaptiveStep(
        &central_icrs_trajectory,
        Ephemeris<ICRS>::NoIntrinsicAcceleration,
        t0 + 2 * *initial_osculating_orbit.elements_at_epoch().period,
        Ephemeris<ICRS>::AdaptiveStepParameters(
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Ephemeris<ICRS>::NewtonianMotionEquation>(),
            /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second)));
    DiscreteTrajectory<GCRS> gcrs_trajectory;
    DiscreteTrajectory<GCRS> ascending;
    DiscreteTrajectory<GCRS> descending;
    for (auto const& [t, dof] : central_icrs_trajectory) {
      CHECK_OK(gcrs_trajectory.Append(t, gcrs_.ToThisFrameAtTime(t)(dof)));
    }
    CHECK_OK(ComputeNodes(gcrs_trajectory,
                          gcrs_trajectory.begin(),
                          gcrs_trajectory.end(),
                          gcrs_trajectory.back().time,
                          Vector<double, GCRS>({0, 0, 1}),
                          /*max_points=*/std::numeric_limits<int>::max(),
                          ascending,
                          descending));
    Instant const t_node =
        ascending
            .upper_bound(
                t0 + *initial_osculating_orbit.elements_at_epoch().period / 2)
            ->time;
    Length δ = 50 * Metre;
    for (int x = -20; x <= 20; ++x) {
      for (int y = -20; y <= 20; ++y) {
        if (x == 0) {
          continue;
        }
        Displacement<LVLH> const r({x * δ, y * δ, 0 * Metre});
        Displacement<LVLH> const r_circular(
            {r.coordinates().x * 2, r.coordinates().y, 0 * Metre});
        if (r_circular.Norm() > 1 * Kilo(Metre)) {
          continue;
        }

        Bivector<double, LVLH> orbit_normal({0, 0, 1});
        auto const circular_tangent = Normalize(orbit_normal * r_circular);
        Vector<double, LVLH> elliptical_tangent(
            {circular_tangent.coordinates().x / 2,
             circular_tangent.coordinates().y,
             circular_tangent.coordinates().z});
        auto const Δa = [&](Speed const v_lvlh) {
          return *KeplerOrbit<ICRS>(
                      earth_,
                      MasslessBody{},
                      lvlh.FromThisFrameAtTime(t0)(
                          {LVLH::origin + r, elliptical_tangent * v_lvlh}) -
                          ephemeris_->trajectory(&earth_)
                              ->EvaluateDegreesOfFreedom(t0),
                      t0)
                      .elements_at_epoch()
                      .semimajor_axis -
                 *initial_osculating_orbit.elements_at_epoch().semimajor_axis;
        };
        absl::btree_set<Speed> v_lvlh_keplerian =
            DoubleBrent(Δa,
                        /*lower_bound=*/-100 * Metre / Second,
                        /*lower_bound=*/100 * Metre / Second);
        CHECK_EQ(v_lvlh_keplerian.size(), 1);

        auto const Δr² = [&](Velocity<LVLH> const& v_lvlh) {
          DiscreteTrajectory<ICRS> trial;
          CHECK_OK(trial.Append(
              t0, lvlh.FromThisFrameAtTime(t0)({LVLH::origin + r, v_lvlh})));
          CHECK_OK(ephemeris_->FlowWithAdaptiveStep(
              &trial,
              Ephemeris<ICRS>::NoIntrinsicAcceleration,
              t_node,
              Ephemeris<ICRS>::AdaptiveStepParameters(
                  EmbeddedExplicitRungeKuttaNyströmIntegrator<
                      DormandالمكاوىPrince1986RKN434FM,
                      Ephemeris<ICRS>::NewtonianMotionEquation>(),
                  /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
                  /*length_integration_tolerance=*/1 * Milli(Metre),
                  /*speed_integration_tolerance=*/1 * Milli(Metre) / Second)));
          return ((lvlh.ToThisFrameAtTime(t_node).rigid_transformation()(
                       trial.EvaluatePosition(t_node)) -
                   LVLH::origin) -
                  r)
              .Norm²();
        };
        auto const grad_Δr² = [&](Velocity<LVLH> const& v_lvlh) {
          Speed const δv = 1 * Micro(Metre) / Second;
          auto const Δr²_v_lvlh = Δr²(v_lvlh);
          return Vector<Quotient<Area, Speed>, LVLH>(
              {(Δr²(v_lvlh +
                    Velocity<LVLH>(
                        {δv, 0 * Metre / Second, 0 * Metre / Second})) -
                Δr²_v_lvlh) /
                   δv,
               (Δr²(v_lvlh +
                    Velocity<LVLH>(
                        {0 * Metre / Second, δv, 0 * Metre / Second})) -
                Δr²_v_lvlh) /
                   δv,
               0 * Metre * Second});
        };

        auto const result = BroydenFletcherGoldfarbShanno<Area, Velocity<LVLH>>(
            /*start_argument=*/elliptical_tangent * *v_lvlh_keplerian.begin(),
            Δr²,
            grad_Δr²,
            10 * Micro(Metre) / Second);
        auto const v_lvlh = result.value();


        LOG(INFO) << "x = " << x << "δ, y = " << y << "δ:\n    v_lvlh = " << v_lvlh
                  << ";\n    Δr² = " << Δr²(v_lvlh);
        icrs_trajectories.emplace_back();
        CHECK_OK(icrs_trajectories.back().Append(
            t0, lvlh.FromThisFrameAtTime(t0)({LVLH::origin + r, v_lvlh})));
        CHECK_OK(ephemeris_->FlowWithAdaptiveStep(
            &icrs_trajectories.back(),
            Ephemeris<ICRS>::NoIntrinsicAcceleration,
            t_node,
            Ephemeris<ICRS>::AdaptiveStepParameters(
                EmbeddedExplicitRungeKuttaNyströmIntegrator<
                    DormandالمكاوىPrince1986RKN434FM,
                    Ephemeris<ICRS>::NewtonianMotionEquation>(),
                /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
                /*length_integration_tolerance=*/1 * Milli(Metre),
                /*speed_integration_tolerance=*/1 * Milli(Metre) / Second)));

        logger.Append("fischbacherField", std::tuple{r, v_lvlh}, ExpressInSIUnits);
      }
    }
    constexpr int n = 50;
    logger.Append("abFormation", std::tuple{0, 0});
    for (auto const& trajectory : icrs_trajectories) {
      DiscreteTrajectory<ICRS> apoapsides;
      DiscreteTrajectory<ICRS> periapsides;
      if (&trajectory == &central_icrs_trajectory) {
        continue;
      }
      ComputeApsides(central_icrs_trajectory,
                     trajectory,
                     trajectory.begin(),
                     trajectory.end(),
                     trajectory.back().time,
                     std::numeric_limits<int>::max(),
                     apoapsides,
                     periapsides);
      logger.Append(
          "abFormation",
          std::tuple{
              (lvlh.ToThisFrameAtTime(apoapsides.front().time)
                   .rigid_transformation()(
                       apoapsides.front().degrees_of_freedom.position()) -
               LVLH::origin)
                  .Norm(),
              (lvlh.ToThisFrameAtTime(periapsides.front().time)
                   .rigid_transformation()(
                       periapsides.front().degrees_of_freedom.position()) -
               LVLH::origin)
                  .Norm()},
          ExpressInSIUnits);
    }
    for (int i = 0; i <= n; ++i) {
      Instant const t = t0 + i * (t_node - t0) / n;
      std::vector<Position<LVLH>> positions;
      for (auto const& trajectory : icrs_trajectories) {
        positions.push_back(lvlh.ToThisFrameAtTime(t).rigid_transformation()(
            trajectory.EvaluatePosition(t)));
      }
      logger.Append("frames", positions, ExpressInSIUnits);
    }
}

TEST_F(FishyTest, Geometric) {
  for (auto const& parameters : std::vector<FishyParameters>{
           {
               .name = "Central",
            .max_degree = 0,
            .zonal_only = true,
            .srp = false,
            .min_t0 = J2000,
              .aspect_ratio = 2,
           },
           {
               .name = "J2",
              .max_degree = 2,
              .zonal_only = true,
              .srp = false,
              .min_t0 = J2000,
              .aspect_ratio = 2,
           },
           {
               .name = "J2AdjustedEllipse",
              .max_degree = 2,
              .zonal_only = true,
              .srp = false,
              .min_t0 = J2000,
              .aspect_ratio = 2 * 1.037,
           },
           {
               .name = "FullGeopotential",
              .max_degree = std::numeric_limits<int>::max(),
              .zonal_only = false,
              .srp = false,
              .min_t0 = J2000,
              .aspect_ratio = 2 / 1.0037,
           },
           {
               .name = "SRP",
              .max_degree = std::numeric_limits<int>::max(),
              .zonal_only = false,
              .srp = true,
              .min_t0 = J2000,
              .aspect_ratio = 2 / 1.0037,
           },
           {
               .name = "FullGeopotentialMay",
              .max_degree = std::numeric_limits<int>::max(),
              .zonal_only = false,
              .srp = false,
              .min_t0 = "2000-05-07T12:00:00"_TT,
              .aspect_ratio = 2 / 1.0037,
           },
           {
               .name = "FullGeopotentialSolstice",
              .max_degree = std::numeric_limits<int>::max(),
              .zonal_only = false,
              .srp = false,
              .min_t0 = "2000-06-21T12:00:00"_TT,
              .aspect_ratio = 2 / 1.0037,
           },
           {
               .name = "SRPMay",
              .max_degree = std::numeric_limits<int>::max(),
              .zonal_only = false,
              .srp = true,
              .min_t0 = "2000-05-07T12:00:00"_TT,
              .aspect_ratio = 2 / 1.0037,
           },
           {
               .name = "SRPSolstice",
              .max_degree = std::numeric_limits<int>::max(),
              .zonal_only = false,
              .srp = true,
              .min_t0 = "2000-06-21T12:00:00"_TT,
              .aspect_ratio = 2 / 1.0037,
           },
       }) {
  SolarSystem<ICRS> solar_system_1950_(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  solar_system_1950_.LimitOblatenessToDegree("Earth", parameters.max_degree);
  if (parameters.zonal_only)
  solar_system_1950_.LimitOblatenessToZonal("Earth");

  not_null<std::unique_ptr<Ephemeris<ICRS>>> ephemeris_ =
      solar_system_1950_.MakeEphemeris(
          /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                   /*geopotential_tolerance=*/0x1p-24},
          Ephemeris<ICRS>::FixedStepParameters(
              SymmetricLinearMultistepIntegrator<
                  QuinlanTremaine1990Order12,
                  Ephemeris<ICRS>::NewtonianMotionEquation>(),
              /*step=*/10 * Minute));
  RotatingBody<ICRS> const& earth_ =
      *solar_system_1950_.rotating_body(*ephemeris_, "Earth");
  RotatingBody<ICRS> const& sun_ =
      *solar_system_1950_.rotating_body(*ephemeris_, "Sun");
  BodyCentredNonRotatingReferenceFrame<ICRS, GCRS> const gcrs_(ephemeris_.get(),
                                                               &earth_);
  for (int month = 2; month <= 12; ++month) {
    LOG(INFO) << "Prolonging to " << std::format("2000-{:02}-01...", month);
    CHECK_OK(ephemeris_->Prolong(
        ParseTT(std::format("2000-{:02}-01T00:00:00", month))));
  }
  Instant const t_max = J2000 + 1 * JulianYear;
  CHECK_OK(ephemeris_->Prolong(t_max));
    double inclination = 97.95; /*
    std::cout << "Enter inclination in degrees ]90, 100[: ";
    std::cin >> inclination;
    if (!(inclination > 90 && inclination < 100)) {
      continue;
    }*/
    double lan = 10;            /*
    std::cout << "Enter longitude of ascending in degrees [0, 360[: ";
    std::cin >> lan;
    if (!(lan >= 0 && lan < 360)) {
      continue;
    }*/
    std::vector<DiscreteTrajectory<ICRS>> icrs_trajectories;
    icrs_trajectories.reserve(400);
    icrs_trajectories.emplace_back();
  auto& central_icrs_trajectory = icrs_trajectories.back();
  KeplerOrbit<ICRS> initial_osculating_orbit(
      earth_,
      MasslessBody{},
      KeplerianElements<ICRS>{
          .eccentricity = 0,
          .semimajor_axis = earth_.mean_radius() + 650 * Kilo(Metre),
          .inclination = inclination * Degree,
          .longitude_of_ascending_node = lan * Degree,
          .argument_of_periapsis = 0 * Degree,
          .mean_anomaly = 0 * Degree,
      },
      J2000);
  CHECK_OK(central_icrs_trajectory.Append(
      J2000,
      ephemeris_->trajectory(&earth_)->EvaluateDegreesOfFreedom(J2000) +
          initial_osculating_orbit.StateVectors(J2000)));
  BodyCentredBodyDirectionReferenceFrame<ICRS, LVLH> lvlh(
      ephemeris_.get(),
      [&]() -> Trajectory<ICRS> const& { return central_icrs_trajectory; },
      &earth_);

  Pressure const P0_srp = TotalSolarIrradiance / SpeedOfLight;
  LOG(INFO) << "P0 = " << P0_srp;
  Mass const m = 575 * Kilogram;
  double const Cr = parameters.srp ? 1.5 : 0;
  Area const A_sat = 105 * Pow<2>(Metre);

  Instant t0 = J2000;
  Length a_osculating_at_t0 =
      *initial_osculating_orbit.elements_at_epoch().semimajor_axis;
  if (parameters.min_t0 != J2000) {
    LOG(INFO) << "Flowing central satellite to " << TTDay(parameters.min_t0)
              << "...";
    CHECK_OK(ephemeris_->FlowWithAdaptiveStep(
        &central_icrs_trajectory,
        [&](Instant const& t,
            DegreesOfFreedom<ICRS> const& dof) -> Vector<Acceleration, ICRS> {
          auto const satellite_sun =
              ephemeris_->trajectory(&sun_)->EvaluatePosition(t) -
              dof.position();
          auto const satellite_earth =
              ephemeris_->trajectory(&earth_)->EvaluatePosition(t) -
              dof.position();
          Length const ray_height =
              satellite_earth.OrthogonalizationAgainst(satellite_sun).Norm() -
              earth_.mean_radius();
          if (ray_height < Length{}) {
            return {};
          } else {
            Pressure const P_srp =
                P0_srp * (Pow<2>(AstronomicalUnit) / satellite_sun.Norm²());
            return -Normalize(satellite_sun) * P_srp * Cr * A_sat / m;
          }
        },
        parameters.min_t0 + 2 * *initial_osculating_orbit.elements_at_epoch().period,
        Ephemeris<ICRS>::GeneralizedAdaptiveStepParameters(
            EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                Fine1987RKNG34,
                Ephemeris<ICRS>::GeneralizedNewtonianMotionEquation>(),
            /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
            /*length_integration_tolerance=*/1 * Centi(Metre),
            /*speed_integration_tolerance=*/1 * Centi(Metre) / Second)));
    DiscreteTrajectory<GCRS> central_gcrs_trajectory;
    for (auto it = central_icrs_trajectory.lower_bound(parameters.min_t0);
         it != central_icrs_trajectory.end();
         ++it) {
      CHECK_OK(central_gcrs_trajectory.Append(
          it->time, gcrs_.ToThisFrameAtTime(it->time)(it->degrees_of_freedom)));
    }
    DiscreteTrajectory<GCRS> ascending;
    DiscreteTrajectory<GCRS> descending;
    CHECK_OK(
        ComputeNodes(central_gcrs_trajectory,
                     central_gcrs_trajectory.upper_bound(parameters.min_t0),
                     central_gcrs_trajectory.end(),
                     central_gcrs_trajectory.back().time,
                     Vector<double, GCRS>({0, 0, 1}),
                     /*max_points=*/std::numeric_limits<int>::max(),
                     ascending,
                     descending));
    t0 = ascending.front().time;
    a_osculating_at_t0 = *KeplerOrbit<ICRS>(
        earth_,
        MasslessBody{},
        central_icrs_trajectory.EvaluateDegreesOfFreedom(t0) -
            ephemeris_->trajectory(&earth_)->EvaluateDegreesOfFreedom(t0),
        t0).elements_at_epoch().semimajor_axis;
  }
  Length δ = 150 * Metre;
  int S1 = -1;
  for (int x = -10; x <= 10; ++x) {
    for (int y = -10; y <= 10; ++y) {
      if (x == 0) {
        continue;
      }
      Displacement<LVLH> const r({x * δ,
                                  y * δ, 0 * Metre});
      Displacement<LVLH> const r_circular(
          {r.coordinates().x * parameters.aspect_ratio, r.coordinates().y, 0 * Metre});
      if (r_circular.Norm() > 1 * Kilo(Metre)) {
        continue;
      }
      Bivector<double, LVLH> orbit_normal({0, 0, 1});
      auto const circular_tangent = Normalize(orbit_normal * r_circular);
      Vector<double, LVLH> elliptical_tangent(
          {circular_tangent.coordinates().x / parameters.aspect_ratio,
           circular_tangent.coordinates().y,
           circular_tangent.coordinates().z});
      auto const Δa = [&](Speed const v_lvlh) {
        return *KeplerOrbit<ICRS>(
                    earth_,
                    MasslessBody{},
                    lvlh.FromThisFrameAtTime(t0)(
                        {LVLH::origin + r, elliptical_tangent * v_lvlh}) -
                        ephemeris_->trajectory(&earth_)
                            ->EvaluateDegreesOfFreedom(t0),
                    t0)
                    .elements_at_epoch()
                    .semimajor_axis -
               a_osculating_at_t0;
      };
      absl::btree_set<Speed> v_lvlh =
          DoubleBrent(Δa,
                      /*lower_bound=*/-100 * Metre / Second,
                      /*lower_bound=*/100 * Metre / Second);
      if (v_lvlh.size() == 0) {
        LOG(INFO) << "x = " << x << "δ, y = " << y
            << "δ: could not find a solution for Δa = 0, with LVLH direction "
            << elliptical_tangent << ". Skipping.";
        continue;
      } else if (v_lvlh.size() > 1) {
        LOG(INFO) << "x = " << x << "δ, y = " << y
                  << "δ: multiple solutions for v_lvlh:";
        for (auto const& v : v_lvlh) {
          LOG(INFO) << "        " << v;
        }
      } else {
        LOG(INFO) << "x = " << x << "δ, y = " << y
                  << "δ: v_lvlh = " << *v_lvlh.begin();
      }
      if (y == 0 && S1 < 0) {
        S1 = icrs_trajectories.size();
      }

      icrs_trajectories.emplace_back();
      CHECK_OK(icrs_trajectories.back().Append(
          t0,
          lvlh.FromThisFrameAtTime(t0)(
              {LVLH::origin + r, elliptical_tangent * *v_lvlh.begin()})));
    }
  }
  LOG(INFO) << "Flowing " << icrs_trajectories.size() << " trajectories...";
  for (int i = 1; i <= 52; ++i) {
    Instant const t = J2000 + i * ((t_max - J2000) / 52);
    if (t < t0) {
      continue;
    }
    Bundle bundle;
    LOG(INFO) << "Flowing to " << TTDay(t);
    for (auto& trajectory : icrs_trajectories) {
      bundle.Add([&]() {
        auto const instance = ephemeris_->NewInstance(
            {&trajectory},
            Ephemeris<ICRS>::NoIntrinsicAccelerations,
            {SymmetricLinearMultistepIntegrator<
                 Quinlan1999Order8A,
                 Ephemeris<ICRS>::NewtonianMotionEquation>(),
             /*step=*/10 * Second});
        return ephemeris_->FlowWithAdaptiveStep(&trajectory, 
            [&](Instant const& t, DegreesOfFreedom<ICRS> const& dof)
                -> Vector<Acceleration, ICRS> {
              auto const satellite_sun =
                  ephemeris_->trajectory(&sun_)->EvaluatePosition(t) -
                  dof.position();
              auto const satellite_earth =
                  ephemeris_->trajectory(&earth_)->EvaluatePosition(t) -
                  dof.position();
              Length const ray_height =
                  satellite_earth.OrthogonalizationAgainst(satellite_sun)
                      .Norm() -
                  earth_.mean_radius();
              if (ray_height < Length{}) {
                return {};
              } else {
                Pressure const P_srp =
                    P0_srp * (Pow<2>(AstronomicalUnit) / satellite_sun.Norm²());
                return -Normalize(satellite_sun) * P_srp * Cr * A_sat / m;
              }
            },
            t,
            Ephemeris<ICRS>::GeneralizedAdaptiveStepParameters(
                EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                    Fine1987RKNG34,
                    Ephemeris<
                        ICRS>::GeneralizedNewtonianMotionEquation>(),
                /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
                /*length_integration_tolerance=*/1 * Centi(Metre),
                /*speed_integration_tolerance=*/1 * Centi(Metre) / Second));
      });
    }
    CHECK_OK(bundle.Join());
  }
  Logger logger(TEMP_DIR / ("fishy" + parameters.name + ".wl"),
                /*make_unique=*/false);
  std::vector <DiscreteTrajectory<GCRS>> gcrs_trajectories;
  std::vector<std::optional<OrbitalElements>> mean_elements;
  std::vector<std::optional<OrbitGroundTrack>> ground_tracks;
  gcrs_trajectories.resize(icrs_trajectories.size());
  mean_elements.resize(icrs_trajectories.size());
  ground_tracks.resize(icrs_trajectories.size());
  Bundle bundle;
  for (int i = 0; i < icrs_trajectories.size(); ++i) {
    bundle.Add([&icrs_trajectories,
                &gcrs_trajectories,
                &ground_tracks,
                &mean_elements,
                &gcrs_,
                &earth_,
                i]() {
    for (auto const& [t, dof] : icrs_trajectories[i]) {
      CHECK_OK(gcrs_trajectories[i].Append(t, gcrs_.ToThisFrameAtTime(t)(dof)));
    }
    PolynomialInMonomialBasis<Angle, Instant, 2> const newcomb_mean_longitude(
        {279 * Degree + 41 * ArcMinute + 48.04 * ArcSecond,
         129'602'768.13 * ArcSecond / (100 * JulianYear),
         1.089 * ArcSecond / Pow<2>(100 * JulianYear)},
        "1899-12-31T12:00:00"_TT);
    OrbitGroundTrack::MeanSun const mean_sun{
        .epoch = J2000,
        .mean_longitude_at_epoch = newcomb_mean_longitude(J2000),
        .year =
            2 * π * Radian / newcomb_mean_longitude.EvaluateDerivative(J2000)};
    ground_tracks[i] =
        OrbitGroundTrack::ForTrajectory(gcrs_trajectories[i], earth_, mean_sun)
            .value();
    mean_elements[i] = OrbitalElements::ForTrajectory(
                           gcrs_trajectories[i], earth_, MasslessBody{})
                           .value();
    LOG(INFO) << "Analysed orbit of #" << i;
      return absl::OkStatus();
                          });
  }
  CHECK_OK(bundle.Join());
  LOG(INFO) << "Sidereal period: " << mean_elements[0]->sidereal_period();
  LOG(INFO) << "Anomalistic period: " << mean_elements[0]->anomalistic_period();
  LOG(INFO) << "Nodal period: " << mean_elements[0]->nodal_period();
  LOG(INFO) << "Initial osculating period: "
            << *initial_osculating_orbit.elements_at_epoch().period;
  LOG(INFO) << "S1 LVLH y at t0: "
            << (lvlh.ToThisFrameAtTime(t0).rigid_transformation()(
                    icrs_trajectories[S1].EvaluatePosition(t0)) -
                LVLH::origin)
                   .coordinates()
                   .y;
  Time const shape_cycle = Brent(
      [&](Time const Δt) {
        return (lvlh.ToThisFrameAtTime(t0 + Δt).rigid_transformation()(
                    icrs_trajectories[S1].EvaluatePosition(t0 + Δt)) -
                LVLH::origin)
            .coordinates()
            .y;
      },
      0.8 * mean_elements[0]->nodal_period(),
      1.2 * mean_elements[0]->nodal_period());
  LOG(INFO) << "Shape cycle: " << shape_cycle;
  Time const sampling_period = shape_cycle;
  std::vector<Instant> times;
  for (int k = 0; k <= 24; ++k) {
    times.push_back(t0 + k * sampling_period / 12);
  }
  for (int n = 3; t0 + n * sampling_period < t_max; ++n) {
    times.push_back(t0 + n * sampling_period);
  }
  logger.Set("times" + parameters.name, times, ExpressInSIUnits);
  for (int i = 0; i < icrs_trajectories.size(); ++i) {
    LOG(INFO) << "Exporting data for trajectory " << i;
    std::vector<Position<LVLH>> lvlh_positions;
    for (Instant const t : times) {
      lvlh_positions.push_back(lvlh.ToThisFrameAtTime(t).rigid_transformation()(icrs_trajectories[i].EvaluatePosition(t)));
    }
    logger.Append(
        "lvlhPositions" + parameters.name, lvlh_positions, ExpressInSIUnits);
    auto const ray_height = [&](Instant const& t) {
      auto const satellite_sun =
          ephemeris_->trajectory(&sun_)->EvaluatePosition(t) -
          icrs_trajectories[i].EvaluatePosition(t);
      auto const satellite_earth =
          ephemeris_->trajectory(&earth_)->EvaluatePosition(t) -
          icrs_trajectories[i].EvaluatePosition(t);
      Length const ray_height =
          satellite_earth.OrthogonalizationAgainst(satellite_sun).Norm() -
          earth_.mean_radius();
      return ray_height;
    };
    std::vector<Instant> eclipse_starts;
    std::vector<Instant> eclipse_ends;
    for (Instant t = t0 + mean_elements[0]->nodal_period() / 2; t < t_max;
         t += mean_elements[0]->nodal_period() / 2) {
      lvlh_positions.push_back(lvlh.ToThisFrameAtTime(t).rigid_transformation()(
          icrs_trajectories[i].EvaluatePosition(t)));
      if (ray_height(t - mean_elements[0]->nodal_period() / 2) < Length{} ||
          ray_height(t) < Length{}) {
        LOG(INFO) << parameters.name <<": Eclipse at node on " << TTDay(t);
        continue;
      }
      Instant const t_min_ray_height = Brent(ray_height, t - mean_elements[0]->nodal_period() / 2, t, std::less<>());
      if (ray_height(t_min_ray_height) < Length{}) {
        eclipse_starts.push_back(Brent(ray_height,
                                       t - mean_elements[0]->nodal_period() / 2,
                                       t_min_ray_height));
        eclipse_ends.push_back(Brent(ray_height, t_min_ray_height, t));
      }
    }
    logger.Append(
        "eclipseStarts" + parameters.name, eclipse_starts, ExpressInSIUnits);
    logger.Append(
        "eclipseEnds" + parameters.name, eclipse_ends, ExpressInSIUnits);
    /*
    logger.Append("meanInclinations",
                  mean_elements[i]->mean_elements() |
                      std::ranges::views::transform(
                          &OrbitalElements::ClassicalElements::inclination) |
                      std::ranges::to<std::vector>(),
                  ExpressInSIUnits);
    logger.Append("meanEccentricities",
                  mean_elements[i]->mean_elements() |
                      std::ranges::views::transform(
                          &OrbitalElements::ClassicalElements::eccentricity) |
                      std::ranges::to<std::vector>(),
                  ExpressInSIUnits);
    logger.Append(
        "meanArgumentsOfPeriapsis",
        mean_elements[i]->mean_elements() |
            std::ranges::views::transform(
                &OrbitalElements::ClassicalElements::argument_of_periapsis) |
            std::ranges::to<std::vector>(),
        ExpressInSIUnits);
    logger.Append("meanElementTimes",
                  mean_elements[i]->mean_elements() |
                      std::ranges::views::transform(
                          &OrbitalElements::ClassicalElements::time) |
                      std::ranges::to<std::vector>(),
                  ExpressInSIUnits);*/
    logger.Append(
        "meanSolarTimesOfAscendingNodes" + parameters.name,
        std::vector{
            ground_tracks[i]->mean_solar_times_of_ascending_nodes()->min,
            ground_tracks[i]->mean_solar_times_of_ascending_nodes()->max},
        ExpressInSIUnits);
    logger.Append(
        "meanSolarTimesOfDescendingNodes" + parameters.name,
        std::vector{
            ground_tracks[i]->mean_solar_times_of_descending_nodes()->min,
            ground_tracks[i]->mean_solar_times_of_descending_nodes()->max},
        ExpressInSIUnits);
  }
  logger.FlushAndClear();
  }
  std::terminate();
}

}  // namespace astronomy
}  // namespace principia