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
#include "numerics/root_finders.hpp"

namespace principia {
namespace astronomy {

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
};

TEST_F(FishyTest, Geometric) {
  for (auto const& parameters : std::vector<FishyParameters>{
           {
               .name = "Central",
            .max_degree = 0,
            .zonal_only = true,
            .srp = false,
           },
           {
               .name = "J2",
              .max_degree = 2,
              .zonal_only = true,
              .srp = false,
           },
           {
               .name = "FullGeopotential",
              .max_degree = std::numeric_limits<int>::max(),
              .zonal_only = false,
              .srp = false,
           },
           {
               .name = "SRP",
              .max_degree = std::numeric_limits<int>::max(),
              .zonal_only = false,
              .srp = true,
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
  CHECK_OK(ephemeris_->Prolong(J2000 + 1 * JulianYear));
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
  Length δ = 150 * Metre;
  int S1 = -1;
  for (int x = -10; x <= 10; ++x) {
    for (int y = -10; y <= 10; ++y) {
      if (x == 0 && y == 0) {
        continue;
      }
      Displacement<LVLH> const r({x * δ,
                                  y * δ, 0 * Metre});
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
                    lvlh.FromThisFrameAtTime(J2000)(
                        {LVLH::origin + r, elliptical_tangent * v_lvlh}) -
                        ephemeris_->trajectory(&earth_)
                            ->EvaluateDegreesOfFreedom(J2000),
                    J2000)
                    .elements_at_epoch()
                    .semimajor_axis -
               *initial_osculating_orbit.elements_at_epoch().semimajor_axis;
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
          J2000,
          lvlh.FromThisFrameAtTime(J2000)(
              {LVLH::origin + r, elliptical_tangent * *v_lvlh.begin()})));
    }
  }
  LOG(INFO) << "Flowing " << icrs_trajectories.size() << " trajectories...";
  for (int i = 1; i <= 52; ++i) {
    Mass const m = 575 * Kilogram;
    double const Cr = parameters.srp ? 1.5 : 0;
    Area const A_sat = 105 * Pow<2>(Metre);
    Instant const t = J2000 + (i / 52.0) * JulianYear;
    Bundle bundle;
    Pressure const P0_srp = TotalSolarIrradiance / SpeedOfLight;
    LOG(INFO) << "P0 = " << P0_srp;
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
  Logger logger(TEMP_DIR / "fishy.wl");
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
  LOG(INFO) << "S1 LVLH y at J2000: "
            << (lvlh.ToThisFrameAtTime(J2000).rigid_transformation()(
                    icrs_trajectories[S1].EvaluatePosition(J2000)) -
                LVLH::origin)
                   .coordinates()
                   .y;
  Time const shape_cycle = Brent(
      [&](Time const Δt) {
        return (lvlh.ToThisFrameAtTime(J2000 + Δt)
                    .rigid_transformation()(
                        icrs_trajectories[S1].EvaluatePosition(J2000 + Δt)) -
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
    times.push_back(J2000 + k * sampling_period / 12);
  }
  for (int n = 3; n * sampling_period < 1 * JulianYear; ++n) {
    times.push_back(J2000 + n * sampling_period);
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
    for (Instant t = J2000 + mean_elements[0]->nodal_period() / 2;
         t < J2000 + 1 * JulianYear;
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