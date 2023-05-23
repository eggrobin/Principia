#include "physics/equipotential.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/plane.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/embedded_explicit_runge_kutta_integrator.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/logger.hpp"
#include "numerics/global_optimization.hpp"
#include "numerics/root_finders.hpp"
#include "physics/body_centred_body_direction_reference_frame.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/ephemeris.hpp"
#include "physics/reference_frame.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace physics {

using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_plane;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::integrators::_embedded_explicit_runge_kutta_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::numerics::_global_optimization;
using namespace principia::numerics::_root_finders;
using namespace principia::physics::_body_centred_body_direction_reference_frame;  // NOLINT
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_equipotential;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_solar_system_factory;

class EquipotentialTest : public ::testing::Test {
 protected:
  using Barycentric = Frame<struct BarycentricTag, Inertial>;
  using World = Frame<struct WorldTag, Arbitrary>;

  EquipotentialTest()
      : ephemeris_parameters_(
            SymmetricLinearMultistepIntegrator<
                QuinlanTremaine1990Order12,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
            /*step=*/10 * Minute),
        solar_system_(make_not_null_unique<SolarSystem<Barycentric>>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt",
            /*ignore_frame=*/true)),
        ephemeris_(solar_system_->MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            ephemeris_parameters_)),
        equipotential_parameters_(
            EmbeddedExplicitRungeKuttaIntegrator<
                DormandPrince1986RK547FC,
                Equipotential<Barycentric, World>::ODE>(),
            /*max_steps=*/1000,
            /*length_integration_tolerance=*/1 * Metre) {}

  Position<World> ComputePositionInWorld(
      Instant const& t,
      ReferenceFrame<Barycentric, World> const& reference_frame,
      SolarSystemFactory::Index const body) {
    auto const to_this_frame = reference_frame.ToThisFrameAtTimeSimilarly(t);
    return to_this_frame.similarity()(
        solar_system_->trajectory(*ephemeris_, SolarSystemFactory::name(body))
            .EvaluatePosition(t));
  }

  std::array<Position<World>, 2> ComputeLagrangePoints(
      SolarSystemFactory::Index const body1,
      SolarSystemFactory::Index const body2,
      Instant const& t,
      ReferenceFrame<Barycentric, World> const& reference_frame,
      Plane<World> const& plane) {
    auto const body1_position =
        ComputePositionInWorld(t, reference_frame, body1);
    auto const body2_position =
        ComputePositionInWorld(t, reference_frame, body2);
    auto const body2_body1 = body1_position - body2_position;

    auto const binormal = plane.UnitBinormals().front();
    Rotation<World, World> const rot_l4(-60 * Degree, binormal);
    auto const body2_l4 = rot_l4(body2_body1);
    auto const l4 = body2_l4 + body2_position;
    Rotation<World, World> const rot_l5(60 * Degree, binormal);
    auto const body2_l5 = rot_l5(body2_body1);
    auto const l5 = body2_l5 + body2_position;

    return {l4, l5};
  }

  // Logs to Mathematica the equipotential line for the given |body| in the
  // specified |reference_frame|.
  void LogEquipotentialLine(
      Logger& logger,
      Plane<World> const& plane,
      Instant const& t,
      ReferenceFrame<Barycentric, World> const& reference_frame,
      SolarSystemFactory::Index const body,
      std::string_view const suffix = "") {
    Equipotential<Barycentric, World> const equipotential(
        equipotential_parameters_,
        &reference_frame,
        /*characteristic_length=*/1 * Metre);
    std::string const name = SolarSystemFactory::name(body);

    CHECK_OK(ephemeris_->Prolong(t));
    auto const line =
        equipotential.ComputeLine(
            plane, t, ComputePositionInWorld(t0_, reference_frame, body));
    std::vector<Position<World>> positions;
    for (auto const& [s, dof] : line) {
      positions.push_back(dof.position());
    }
    logger.Set(absl::StrCat("equipotential", name, suffix),
               positions,
               ExpressIn(Metre));
  }

  // Logs to Mathematica a family of equipotential lines determined by a
  // parameter.  There must exist an overload of |ComputeLine| with a
  // |LineParameter| as its third argument.
  template<typename LineParameter>
  void LogFamilyOfEquipotentialLines(
      Logger& logger,
      ReferenceFrame<Barycentric, World> const& reference_frame,
      int const number_of_days,
      std::string_view const suffix,
      std::function<std::vector<LineParameter>(
          Position<World> const& l4,
          Position<World> const& l5)> const& get_line_parameters) {
    Equipotential<Barycentric, World> const equipotential(
        equipotential_parameters_,
        &reference_frame,
        /*characteristic_length=*/1 * Metre);
    auto const plane =
        Plane<World>::OrthogonalTo(Vector<double, World>({0, 0, 1}));

    std::vector<std::vector<std::vector<Position<World>>>> all_positions;
    for (int j = 0; j < number_of_days; ++j) {
      Instant const t = t0_ + j * Day;
      CHECK_OK(ephemeris_->Prolong(t));
      all_positions.emplace_back();

      auto const& [l4, l5] = ComputeLagrangePoints(SolarSystemFactory::Earth,
                                                    SolarSystemFactory::Moon,
                                                    t,
                                                    reference_frame,
                                                    plane);

      for (auto const& line_parameter : get_line_parameters(l4, l5)) {
        auto const line =
            equipotential.ComputeLine(plane, t, line_parameter);
        all_positions.back().emplace_back();
        for (auto const& [s, dof] : line) {
          all_positions.back().back().push_back(dof.position());
        }
      }
    }
    logger.Set(absl::StrCat("equipotentialsEarthMoon", suffix),
               all_positions,
               ExpressIn(Metre));
  }

  Instant const t0_;
  Ephemeris<Barycentric>::FixedStepParameters const ephemeris_parameters_;
  not_null<std::unique_ptr<SolarSystem<Barycentric>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
  Equipotential<Barycentric, World>::AdaptiveParameters const
      equipotential_parameters_;
};

#if !_DEBUG
TEST_F(EquipotentialTest, BodyCentredNonRotating) {
  Logger logger(TEMP_DIR / "equipotential_bcnr.wl",
                /*make_unique=*/false);
  auto const reference_frame(
      BodyCentredNonRotatingReferenceFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Sun))));
  Equipotential<Barycentric, World> const equipotential(
      equipotential_parameters_,
      &reference_frame,
      /*characteristic_length=*/1 * Metre);

  auto const plane =
      Plane<World>::OrthogonalTo(Vector<double, World>({2, 3, -5}));

  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       reference_frame,
                       SolarSystemFactory::Mercury);
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       reference_frame,
                       SolarSystemFactory::Earth);
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 1 * Day,
                       reference_frame,
                       SolarSystemFactory::Jupiter, "Close");
  LogEquipotentialLine(logger,
                       plane,
                       t0_ + 100 * Day,
                       reference_frame,
                       SolarSystemFactory::Jupiter, "Far");
}

TEST_F(EquipotentialTest, BodyCentredBodyDirection_EquidistantPoints) {
  Logger logger(TEMP_DIR / "equipotential_bcbd_distances.wl",
                /*make_unique=*/false);
  auto const reference_frame(
      BodyCentredBodyDirectionReferenceFrame<Barycentric, World>(
          ephemeris_.get(),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Earth)),
          solar_system_->massive_body(
              *ephemeris_,
              SolarSystemFactory::name(SolarSystemFactory::Moon))));

  LogFamilyOfEquipotentialLines<Position<World>>(
      logger,
      reference_frame,
      /*number_of_days=*/30,
      /*suffix=*/"Distances",
      [](Position<World> const& l4, Position<World> const& l5) {
        std::vector<Position<World>> positions;
        for (int i = 0; i <= 10; ++i) {
          positions.push_back(Barycentre(
              std::pair{l4, l5}, std::pair{i / 10.0, (10.0 - i) / 10.0}));
        }
        return positions;
      });
}

TEST_F(EquipotentialTest, DISABLED_RotatingPulsating_GlobalOptimization) {
  Logger logger(TEMP_DIR / "equipotential_rp_global.wl",
                /*make_unique=*/false);
  std::int64_t const number_of_days = 502;
  auto const earth = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Earth));
  auto const moon = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Moon));
  auto const reference_frame(
      RotatingPulsatingReferenceFrame<Barycentric, World>(
          ephemeris_.get(), moon, earth));
  CHECK_OK(ephemeris_->Prolong(t0_ + number_of_days * Day));

  DegreesOfFreedom<Barycentric> const earth_dof =
      ephemeris_->trajectory(earth)->EvaluateDegreesOfFreedom(t0_);
  DegreesOfFreedom<Barycentric> const moon_dof =
      ephemeris_->trajectory(moon)->EvaluateDegreesOfFreedom(t0_);
  KeplerOrbit<Barycentric> const moon_orbit(
      *earth,
      *moon,
      moon_dof - earth_dof,
      t0_);
  KeplerianElements<Barycentric> const moon_elements =
      moon_orbit.elements_at_epoch();

  KeplerianElements<Barycentric> const elements{
      .periapsis_distance = 71'000 * Kilo(Metre),
      .apoapsis_distance = 0.65 * moon_elements.periapsis_distance.value(),
      .inclination = moon_elements.inclination,
      .longitude_of_ascending_node = moon_elements.longitude_of_ascending_node,
      .argument_of_periapsis = *moon_elements.argument_of_periapsis + Degree,
      .mean_anomaly = 0 * Degree};
  auto const earth_world_dof =
      reference_frame.ToThisFrameAtTimeSimilarly(t0_)(earth_dof);
  auto const moon_world_dof =
      reference_frame.ToThisFrameAtTimeSimilarly(t0_)(moon_dof);
  Position<World> const q_earth = earth_world_dof.position();
  Position<World> const q_moon = moon_world_dof.position();
  Velocity<World> const v_earth = earth_world_dof.velocity();
  Velocity<World> const v_moon = moon_world_dof.velocity();
  Position<World> const initial_earth_moon_l5 =
      Barycentre(std::pair(q_earth, q_moon), std::pair(1, 1)) +
      (q_earth - q_moon).Norm() *
          Vector<double, World>({0, quantities::Sqrt(3) / 2, 0});
  using MEO = Frame<struct MEOTag, Arbitrary>;
  BodyCentredBodyDirectionReferenceFrame<Barycentric, MEO> meo(
      ephemeris_.get(), moon, earth);
  // The initial states for four trajectories:
  // [0]: initially stationary in the rotating-pulsating frame near L3;
  // [1]: initially stationary in MEO at L5;
  // [2]: initially stationary in the rotating-pulsating frame at L5;
  // [3]: in an elliptic Earth orbit that reaches 65% of the way to the Moon.
  std::vector<DegreesOfFreedom<Barycentric>> const initial_states{
      reference_frame.FromThisFrameAtTimeSimilarly(t0_)(
          {q_earth + (q_earth - q_moon), World::unmoving}),
      meo.FromThisFrameAtTime(t0_)(
          {meo.ToThisFrameAtTime(t0_).rigid_transformation()(
               reference_frame.FromThisFrameAtTimeSimilarly(t0_).similarity()(
                   initial_earth_moon_l5)),
           MEO::unmoving}),
      reference_frame.FromThisFrameAtTimeSimilarly(t0_)(
          {initial_earth_moon_l5, World::unmoving}),
      earth_dof +
          KeplerOrbit<Barycentric>(*earth, MasslessBody{}, elements, t0_)
              .StateVectors(t0_)};

  std::vector<not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>>>
      trajectories;
  std::vector<not_null<DiscreteTrajectory<Barycentric>*>> instance_trajectories;
  for (auto const& s : initial_states) {
    trajectories.push_back(
        make_not_null_unique<DiscreteTrajectory<Barycentric>>());
    instance_trajectories.push_back(trajectories.back().get());
    CHECK_OK(trajectories.back()->Append(t0_, s));
    trajectories.back()->segments().front().SetDownsampling(
        {.max_dense_intervals = 10'000, .tolerance = 10 * Metre});
  }
  auto const instance = ephemeris_->NewInstance(
      instance_trajectories,
      Ephemeris<Barycentric>::NoIntrinsicAccelerations,
      Ephemeris<Barycentric>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              Quinlan1999Order8A,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*step=*/10 * Second));

  LOG(ERROR) << "Flowing trajectories";
  CHECK_OK(
      ephemeris_->FlowWithFixedStep(t0_ + number_of_days * Day, *instance));
  LOG(ERROR) << "Flowed";

  Instant t = t0_;
  auto const potential = [&reference_frame,
                          &t](Position<World> const& position) {
    return reference_frame.GeometricPotential(t, position);
  };
  auto const acceleration = [&reference_frame,
                              &t](Position<World> const& position) {
    auto const acceleration =
        reference_frame.GeometricAcceleration(t, {position, Velocity<World>{}});
    // Note the sign.
    return -Vector<Acceleration, World>({acceleration.coordinates()[0],
                                         acceleration.coordinates()[1],
                                         Acceleration{}});
  };
  const MultiLevelSingleLinkage<SpecificEnergy, Position<World>, 2>::Box box = {
      .centre = World::origin,
      .vertices = {Displacement<World>({3 * Metre,
                                        0 * Metre,
                                        0 * Metre}),
                   Displacement<World>({0 * Metre,
                                        3 * Metre,
                                        0 * Metre})}};

  MultiLevelSingleLinkage<SpecificEnergy, Position<World>, 2> optimizer(
      box, potential, acceleration);
  constexpr Length characteristic_length = 1 * Nano(Metre);
  Equipotential<Barycentric, World> const equipotential(
      {EmbeddedExplicitRungeKuttaIntegrator<
           DormandPrince1986RK547FC,
           Equipotential<Barycentric, World>::ODE>(),
       /*max_steps=*/1000,
       /*length_integration_tolerance=*/characteristic_length},
      &reference_frame,
      characteristic_length);
  auto const plane =
      Plane<World>::OrthogonalTo(Vector<double, World>({0, 0, 1}));

  std::vector<std::vector<std::vector<std::vector<Position<World>>>>>
      all_positions;
  std::vector<std::vector<Position<World>>> trajectory_positions(
      trajectories.size());
  for (int j = 0; j < number_of_days; ++j) {
    LOG(ERROR) << "Day #" << j;
    t = t0_ + j * Day;
    CHECK_OK(ephemeris_->Prolong(t));
    std::vector<std::vector<std::vector<Position<World>>>>&
        equipotentials_at_t = all_positions.emplace_back();

    auto const arg_maximorum = optimizer.FindGlobalMaxima(
        /*points_per_round=*/1000,
        /*number_of_rounds=*/std::nullopt,
        /*local_search_tolerance=*/1e-3 * Metre);
    logger.Append("argMaximorum", arg_maximorum, ExpressIn(Metre));
    std::vector<SpecificEnergy> maxima;
    SpecificEnergy maximum_maximorum = -Infinity<SpecificEnergy>;

    Position<World> const earth_position =
        reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
            ephemeris_->trajectory(earth)->EvaluatePosition(t));
    Position<World> const moon_position =
        reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
            ephemeris_->trajectory(moon)->EvaluatePosition(t));
    for (auto const& arg_maximum : arg_maximorum) {
      maxima.push_back(potential(arg_maximum));
      maximum_maximorum = std::max(maximum_maximorum, maxima.back());
    }
    logger.Append("maxima", maxima, ExpressIn(Metre, Second));

    double const arg_approx_l1 = Brent(
        [&](double const x) {
          return potential(Barycentre(std::pair(moon_position, earth_position),
                                      std::pair(x, 1 - x)));
        },
        0.0,
        1.0,
        std::greater<>{});
    double const arg_approx_l2 = Brent(
        [&](double x) {
          return potential(Barycentre(
              std::pair(moon_position,
                        World::origin + 2 * (moon_position - World::origin)),
              std::pair(x, 1 - x)));
        },
        0.0,
        1.0,
        std::greater<>{});
    SpecificEnergy const approx_l1_energy =
        potential(Barycentre(std::pair(moon_position, earth_position),
                             std::pair(arg_approx_l1, 1 - arg_approx_l1)));
    SpecificEnergy const approx_l2_energy = potential(Barycentre(
        std::pair(moon_position,
                  World::origin + 2 * (moon_position - World::origin)),
        std::pair(arg_approx_l2, 1 - arg_approx_l2)));
    for (int i = 0; i < trajectories.size(); ++i) {
      DegreesOfFreedom<World> const dof =
          reference_frame.ToThisFrameAtTimeSimilarly(t)(
              trajectories[i]->EvaluateDegreesOfFreedom(t));
      trajectory_positions[i].push_back(dof.position());
    }
    // TODO(egg): Somehow extract that from the reference frame.
    auto const r = [&](Instant const& t) -> Length {
      return (ephemeris_->trajectory(earth)->EvaluatePosition(t) -
              ephemeris_->trajectory(moon)->EvaluatePosition(t))
          .Norm();
    };
    SpecificEnergy const ΔV = maximum_maximorum - approx_l1_energy;
    for (int i = 1; i <= 8; ++i) {
      SpecificEnergy const energy = maximum_maximorum - i * (1.0 / 7.0 * ΔV);
      std::vector<std::vector<Position<World>>>& equipotentials_at_energy =
          equipotentials_at_t.emplace_back();
      auto const& lines = equipotential.ComputeLines(
          plane,
          t,
          arg_maximorum,
          {{moon_position, moon->min_radius() / r(t) * (1 * Metre)},
           {earth_position, earth->min_radius() / r(t) * (1 * Metre)}},
          [](Position<World> q) {
            return World::origin + Normalize(q - World::origin) * 3 * Metre;
          },
          energy);
      for (auto const& line : lines) {
        std::vector<Position<World>>& equipotential =
            equipotentials_at_energy.emplace_back();
        for (auto const& [s, dof] : line) {
          equipotential.push_back(dof.position());
        }
      }
    }
    auto const& lines = equipotential.ComputeLines(
        plane,
        t,
        arg_maximorum,
        {{moon_position, moon->min_radius() / r(t) * (1 * Metre)},
         {earth_position, earth->min_radius() / r(t) * (1 * Metre)}},
        [](Position<World> q) {
          return World::origin + Normalize(q - World::origin) * 3 * Metre;
        },
        approx_l2_energy);
    std::vector<std::vector<Position<World>>>& equipotentials_at_l2_energy =
        equipotentials_at_t.emplace_back();
    for (auto const& line : lines) {
      std::vector<Position<World>>& equipotential =
          equipotentials_at_l2_energy.emplace_back();
      for (auto const& [s, dof] : line) {
        equipotential.push_back(dof.position());
      }
    }
  }
  std::vector<std::vector<Position<World>>> world_trajectories;
  for (auto const& trajectory : trajectories) {
    world_trajectories.emplace_back();
    for (auto const& [t, dof] : *trajectory) {
      world_trajectories.back().push_back(
          reference_frame.ToThisFrameAtTimeSimilarly(t).similarity()(
              dof.position()));
    }
  }
  logger.Set("trajectories", world_trajectories, ExpressIn(Metre));
  logger.Set("trajectoryPositions",
             trajectory_positions,
             ExpressIn(Metre));
  logger.Set("equipotentialsEarthMoonGlobalOptimization",
             all_positions,
             ExpressIn(Metre));
}

TEST_F(EquipotentialTest, RotatingPulsating_SunNeptune) {
  Logger logger(TEMP_DIR / "equipotential_rp_global.wl",
                /*make_unique=*/false);
  auto const sun = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Sun));
  auto const mercury = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Mercury));
  auto const venus = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Venus));
  auto const earth = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Earth));
  auto const moon = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Moon));
  auto const mars = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Mars));
  auto const jupiter = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Jupiter));
  auto const saturn = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Saturn));
  auto const uranus = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Uranus));
  auto const miranda = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Miranda));
  auto const ariel = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Ariel));
  auto const umbriel = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Umbriel));
  auto const titania = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Titania));
  auto const oberon = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Oberon));
  auto const neptune = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Neptune));
  auto const triton = solar_system_->massive_body(
      *ephemeris_, SolarSystemFactory::name(SolarSystemFactory::Triton));
  std::vector<not_null<MassiveBody const*>> all_inner;
  for (int i = 0; i <= SolarSystemFactory::LastBody; ++i) {
    switch (i) {
      case SolarSystemFactory::Uranus:
      case SolarSystemFactory::Miranda:
      case SolarSystemFactory::Ariel:
      case SolarSystemFactory::Umbriel:
      case SolarSystemFactory::Titania:
      case SolarSystemFactory::Oberon:
      case SolarSystemFactory::Neptune:
      case SolarSystemFactory::Triton:
      case SolarSystemFactory::Pluto:
      case SolarSystemFactory::Charon:
      case SolarSystemFactory::Eris:
        continue;
      default:
        all_inner.push_back(solar_system_->massive_body(
            *ephemeris_, SolarSystemFactory::name(i)));
    }
  }
  std::vector<std::vector<std::vector<std::vector<Position<World>>>>>
      equipotentials;
  std::vector<AngularVelocity<Barycentric>> angular_velocities;
  std::vector<AngularVelocity<World>> angular_velocities_in_world;
  std::vector<std::vector<SpecificEnergy>> energies;
  for (auto const [primaries, secondaries] :
       {std::pair{std::vector{sun}, std::vector{uranus}},
        std::pair{
            std::vector{sun},
            std::vector{uranus, miranda, ariel, umbriel, titania, oberon}},
        /*std::pair{
            std::vector{
                sun, mercury, venus, earth, moon, mars, jupiter, saturn},
            std::vector{uranus, miranda, ariel, umbriel, titania, oberon}},*/
        std::pair{
            all_inner,
            std::vector{uranus, miranda, ariel, umbriel, titania, oberon}}}) {
    energies.emplace_back();
    auto const reference_frame(
        RotatingPulsatingReferenceFrame<Barycentric, World>(
            ephemeris_.get(), primaries, secondaries));
    CHECK_OK(ephemeris_->Prolong(t0_));
    angular_velocities.push_back(reference_frame.ToThisFrameAtTimeSimilarly(t0_)
                                     .angular_velocity_of<World>());
    angular_velocities_in_world.push_back(
        reference_frame.ToThisFrameAtTimeSimilarly(t0_)
            .angular_velocity_of<Barycentric>());

    auto const potential = [this,
                            &reference_frame](Position<World> const& position) {
      return reference_frame.GeometricPotential(t0_, position);
    };
    auto const acceleration =
        [this, &reference_frame](Position<World> const& position) {
          auto const acceleration = reference_frame.GeometricAcceleration(
              t0_, {position, Velocity<World>{}});
          // Note the sign.
          return -Vector<Acceleration, World>({acceleration.coordinates()[0],
                                               acceleration.coordinates()[1],
                                               Acceleration{}});
        };
    const MultiLevelSingleLinkage<SpecificEnergy, Position<World>, 2>::Box box =
        {.centre = World::origin,
         .vertices = {Displacement<World>({3 * Metre, 0 * Metre, 0 * Metre}),
                      Displacement<World>({0 * Metre, 3 * Metre, 0 * Metre})}};

    MultiLevelSingleLinkage<SpecificEnergy, Position<World>, 2> optimizer(
        box, potential, acceleration);
    constexpr Length characteristic_length = 1 * Nano(Metre);
    Equipotential<Barycentric, World> const equipotential(
        {EmbeddedExplicitRungeKuttaIntegrator<
             DormandPrince1986RK547FC,
             Equipotential<Barycentric, World>::ODE>(),
         /*max_steps=*/1000,
         /*length_integration_tolerance=*/characteristic_length},
        &reference_frame,
        characteristic_length);
    auto const plane =
        Plane<World>::OrthogonalTo(Vector<double, World>({0, 0, 1}));

    CHECK_OK(ephemeris_->Prolong(t0_));
    std::vector<std::vector<std::vector<Position<World>>>>&
        equipotentials_at_t = equipotentials.emplace_back();

    auto const arg_maximorum = optimizer.FindGlobalMaxima(
        /*points_per_round=*/1000,
        /*number_of_rounds=*/std::nullopt,
        /*local_search_tolerance=*/1e-3 * Metre);
    logger.Append("argMaximorum", arg_maximorum, ExpressIn(Metre));
    std::vector<SpecificEnergy> maxima;
    SpecificEnergy maximum_maximorum = -Infinity<SpecificEnergy>;
    SpecificEnergy minimum_maximorum = Infinity<SpecificEnergy>;

    Position<World> const earth_position =
        reference_frame.ToThisFrameAtTimeSimilarly(t0_).similarity()(
            ephemeris_->trajectory(sun)->EvaluatePosition(t0_));
    Position<World> const moon_position =
        reference_frame.ToThisFrameAtTimeSimilarly(t0_).similarity()(
            ephemeris_->trajectory(uranus)->EvaluatePosition(t0_));
    for (auto const& arg_maximum : arg_maximorum) {
      maxima.push_back(potential(arg_maximum));
      maximum_maximorum = std::max(maximum_maximorum, maxima.back());
      minimum_maximorum = std::min(minimum_maximorum, maxima.back());
    }
    logger.Append("maxima", maxima, ExpressIn(Metre, Second));

    double const arg_approx_l1 = Brent(
        [&](double const x) {
          return potential(Barycentre(std::pair(moon_position, earth_position),
                                      std::pair(x, 1 - x)));
        },
        0.0,
        1.0,
        std::greater<>{});
    double const arg_approx_l2 = Brent(
        [&](double x) {
          return potential(Barycentre(
              std::pair(moon_position,
                        World::origin + 2 * (moon_position - World::origin)),
              std::pair(x, 1 - x)));
        },
        0.0,
        1.0,
        std::greater<>{});
    SpecificEnergy const approx_l1_energy =
        potential(Barycentre(std::pair(moon_position, earth_position),
                             std::pair(arg_approx_l1, 1 - arg_approx_l1)));
    SpecificEnergy const approx_l2_energy = potential(Barycentre(
        std::pair(moon_position,
                  World::origin + 2 * (moon_position - World::origin)),
        std::pair(arg_approx_l2, 1 - arg_approx_l2)));
    logger.Append("approxL1Energy", approx_l1_energy, ExpressIn(Metre, Second));
    logger.Append("approxL2Energy", approx_l2_energy, ExpressIn(Metre, Second));
    // TODO(egg): Somehow extract that from the reference frame.
    auto const r = [&](Instant const& t) -> Length {
      return (ephemeris_->trajectory(sun)->EvaluatePosition(t) -
              ephemeris_->trajectory(uranus)->EvaluatePosition(t))
          .Norm();
    };
    SpecificEnergy ΔV =
        (maximum_maximorum - approx_l1_energy) /
        (4 *
         Sqrt(
             reference_frame.primaries().front()->gravitational_parameter() /
             reference_frame.secondaries().front()->gravitational_parameter()));
    SpecificEnergy const ΔV_max = (maximum_maximorum - approx_l1_energy) / 7;
    SpecificEnergy energy = maximum_maximorum;
    double β = e;
    for (int i = 0; energy >= approx_l1_energy; ++i) {
      if (ΔV * e <= ΔV_max) {
        ΔV *= e;
      }
      energy -= ΔV;
      energies.back().push_back(energy);
      std::vector<std::vector<Position<World>>>& equipotentials_at_energy =
          equipotentials_at_t.emplace_back();
      auto const& lines = equipotential.ComputeLines(
          plane,
          t0_,
          arg_maximorum,
          {{moon_position, uranus->min_radius() / r(t0_) * (1 * Metre)},
           {earth_position, sun->min_radius() / r(t0_) * (1 * Metre)}},
          [](Position<World> q) {
            return World::origin + Normalize(q - World::origin) * 3 * Metre;
          },
          energy);
      for (auto const& line : lines) {
        std::vector<Position<World>>& equipotential =
            equipotentials_at_energy.emplace_back();
        for (auto const& [s, dof] : line) {
          equipotential.push_back(dof.position());
        }
      }
    }
    for (SpecificEnergy const energy : {approx_l1_energy, approx_l2_energy}) {
      auto& equipotentials_at_energy = equipotentials_at_t.emplace_back();
      auto lines = equipotential.ComputeLines(
          plane,
          t0_,
          arg_maximorum,
          {{moon_position, uranus->min_radius() / r(t0_) * (1 * Metre)},
           {earth_position, sun->min_radius() / r(t0_) * (1 * Metre)}},
          [](Position<World> q) {
            return World::origin + Normalize(q - World::origin) * 3 * Metre;
          },
          energy);
      for (auto const& line : lines) {
        std::vector<Position<World>>& equipotential =
            equipotentials_at_energy.emplace_back();
        for (auto const& [s, dof] : line) {
          equipotential.push_back(dof.position());
        }
      }
    }
  }
  logger.Set("energies", energies, ExpressIn(Metre, Second));
  logger.Set("angularVelocitiesSunNeptune",
             angular_velocities,
             ExpressIn(Radian, Second));
  logger.Set("angularVelocitiesInWorldSunNeptune",
             angular_velocities_in_world,
             ExpressIn(Radian, Second));
  logger.Set("equipotentialsSunNeptune", equipotentials, ExpressIn(Metre));
}

#endif

}  // namespace physics
}  // namespace principia
