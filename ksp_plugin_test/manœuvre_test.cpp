﻿
#include <cmath>
#include <functional>
#include <sstream>
#include <type_traits>

#include "geometry/affine_map.hpp"
#include "geometry/frame.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/point.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock-actions.h"
#include "gmock/gmock-generated-matchers.h"
#include "gmock/gmock-generated-nice-strict.h"
#include "gmock/gmock-matchers.h"
#include "gmock/gmock-spec-builders.h"
#include "gtest/gtest.h"
#include "ksp_plugin/manœuvre.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/mock_dynamic_frame.hpp"
#include "physics/rigid_motion.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/si.hpp"
#include "quantities/uk.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using geometry::Frame;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::MockDynamicFrame;
using physics::RigidTransformation;
using quantities::Pow;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using quantities::uk::Foot;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::RelativeError;
using ::testing::AllOf;
using ::testing::Gt;
using ::testing::Lt;
using ::testing::Return;
using ::testing::StrictMock;

namespace ksp_plugin {

class ManœuvreTest : public ::testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;
  using Rendering = Frame<serialization::Frame::TestTag,
                          serialization::Frame::TEST2, false>;

  not_null<std::unique_ptr<StrictMock<MockDynamicFrame<World, Rendering>>>>
  MakeMockDynamicFrame() {
    auto owned_mock_dynamic_frame =
        make_not_null_unique<StrictMock<MockDynamicFrame<World, Rendering>>>();
    mock_dynamic_frame_ = owned_mock_dynamic_frame.get();
    return std::move(owned_mock_dynamic_frame);
  }

  StrictMock<MockDynamicFrame<World, Rendering>> const* mock_dynamic_frame_;
  DiscreteTrajectory<World> discrete_trajectory_;
  DegreesOfFreedom<World> const dof_ = {
      World::origin + Displacement<World>({1 * Metre, 9 * Metre, 5 * Metre}),
      Velocity<World>(
          {8 * Metre / Second, 10 * Metre / Second, 4 * Metre / Second})};
  DegreesOfFreedom<Rendering> const rendering_dof_ = {
      Rendering::origin +
          Displacement<Rendering>({1 * Metre, 9 * Metre, 5 * Metre}),
      Velocity<Rendering>(
          {8 * Metre / Second, 10 * Metre / Second, 4 * Metre / Second})};
  RigidMotion<World, Rendering> const rigid_motion_ =
      RigidMotion<World, Rendering>(
          RigidTransformation<World, Rendering>(
              World::origin,
              Rendering::origin,
              OrthogonalMap<World, Rendering>::Identity()),
          AngularVelocity<World>(), Velocity<World>());
};

TEST_F(ManœuvreTest, TimedBurn) {
  Instant const t0 = Instant();
  Vector<double, Frenet<Rendering>> e_y({0, 1, 0});

  Manœuvre<World, Rendering> manœuvre(
      1 * Newton /*thrust*/,
      2 * Kilogram /*initial_mass*/,
      1 * Newton * Second / Kilogram /*specific_impulse*/,
      2 * e_y /*direction*/,
      MakeMockDynamicFrame());
  EXPECT_EQ(1 * Newton, manœuvre.thrust());
  EXPECT_EQ(2 * Kilogram, manœuvre.initial_mass());
  EXPECT_EQ(1 * Metre / Second, manœuvre.specific_impulse());
  EXPECT_EQ(e_y, manœuvre.direction());
  EXPECT_EQ(1 * Kilogram / Second, manœuvre.mass_flow());

  manœuvre.set_duration(1 * Second);
  EXPECT_EQ(1 * Second, manœuvre.duration());
  EXPECT_EQ(std::log(2) * Metre / Second, manœuvre.Δv());
  EXPECT_EQ((2 - Sqrt(2)) * Second, manœuvre.time_to_half_Δv());
  EXPECT_EQ(1 * Kilogram, manœuvre.final_mass());

  manœuvre.set_initial_time(t0);
  EXPECT_EQ(t0, manœuvre.initial_time());
  EXPECT_EQ(t0 + 1 * Second, manœuvre.final_time());
  EXPECT_EQ(t0 + (2 - Sqrt(2)) * Second, manœuvre.time_of_half_Δv());

  discrete_trajectory_.Append(manœuvre.initial_time(), dof_);
  EXPECT_CALL(*mock_dynamic_frame_, ToThisFrameAtTime(manœuvre.initial_time()))
      .WillOnce(Return(rigid_motion_));
  EXPECT_CALL(*mock_dynamic_frame_,
              FrenetFrame(manœuvre.initial_time(), rendering_dof_))
      .WillOnce(
          Return(Rotation<Frenet<Rendering>, Rendering>::Identity()));
  auto const acceleration = manœuvre.acceleration(discrete_trajectory_);
  EXPECT_EQ(
      0 * Metre / Pow<2>(Second),
      acceleration(manœuvre.initial_time() - 1 * Second).Norm());
  EXPECT_EQ(0.5 * Metre / Pow<2>(Second),
            acceleration(manœuvre.initial_time()).Norm());
  EXPECT_THAT(acceleration(manœuvre.time_of_half_Δv()).Norm(),
              AlmostEquals(Sqrt(0.5) * Metre / Pow<2>(Second), 1));
  EXPECT_EQ(1 * Metre / Pow<2>(Second),
            acceleration(manœuvre.final_time()).Norm());
  EXPECT_THAT(
      acceleration(manœuvre.final_time()),
      Componentwise(0 * Metre / Pow<2>(Second), 1 * Metre / Pow<2>(Second),
                    0 * Metre / Pow<2>(Second)));
  EXPECT_EQ(0 * Metre / Pow<2>(Second),
            acceleration(manœuvre.final_time() + 1 * Second).Norm());
}

TEST_F(ManœuvreTest, TargetΔv) {
  Instant const t0 = Instant();
  Vector<double, Frenet<Rendering>> e_y({0, 1, 0});
  Manœuvre<World, Rendering> manœuvre(
      1 * Newton /*thrust*/,
      2 * Kilogram /*initial_mass*/,
      1 * Newton * Second / Kilogram /*specific_impulse*/,
      e_y /*direction*/,
      MakeMockDynamicFrame());
  EXPECT_EQ(1 * Newton, manœuvre.thrust());
  EXPECT_EQ(2 * Kilogram, manœuvre.initial_mass());
  EXPECT_EQ(1 * Metre / Second, manœuvre.specific_impulse());
  EXPECT_EQ(e_y, manœuvre.direction());
  EXPECT_EQ(1 * Kilogram / Second, manœuvre.mass_flow());

  manœuvre.set_Δv(1 * Metre / Second);
  EXPECT_EQ(1 * Metre / Second, manœuvre.Δv());
  EXPECT_EQ((2 - 2 / e) * Second, manœuvre.duration());
  EXPECT_EQ((2 - 2 / Sqrt(e)) * Second, manœuvre.time_to_half_Δv());
  EXPECT_EQ((2 / e) * Kilogram, manœuvre.final_mass());

  manœuvre.set_time_of_half_Δv(t0);
  EXPECT_EQ(t0 - (2 - 2 / Sqrt(e)) * Second, manœuvre.initial_time());
  EXPECT_EQ(t0 + (2 / Sqrt(e) - 2 / e) * Second, manœuvre.final_time());
  EXPECT_EQ(t0, manœuvre.time_of_half_Δv());

  discrete_trajectory_.Append(manœuvre.initial_time(), dof_);
  EXPECT_CALL(*mock_dynamic_frame_, ToThisFrameAtTime(manœuvre.initial_time()))
      .WillOnce(Return(rigid_motion_));
  EXPECT_CALL(*mock_dynamic_frame_,
              FrenetFrame(manœuvre.initial_time(), rendering_dof_))
      .WillOnce(
          Return(Rotation<Frenet<Rendering>, Rendering>::Identity()));
  auto const acceleration = manœuvre.acceleration(discrete_trajectory_);
  EXPECT_EQ(
      0 * Metre / Pow<2>(Second),
      acceleration(manœuvre.initial_time() - 1 * Second).Norm());
  EXPECT_EQ(0.5 * Metre / Pow<2>(Second),
            acceleration(manœuvre.initial_time()).Norm());
  EXPECT_EQ(Sqrt(e) / 2 * Metre / Pow<2>(Second),
            acceleration(manœuvre.time_of_half_Δv()).Norm());
  EXPECT_EQ((e / 2) * Metre / Pow<2>(Second),
            acceleration(manœuvre.final_time()).Norm());
  EXPECT_EQ(0 * Metre / Pow<2>(Second),
            acceleration(manœuvre.final_time() + 1 * Second).Norm());
}

TEST_F(ManœuvreTest, Apollo8SIVB) {
  // Data from NASA's Saturn V Launch Vehicle, Flight Evaluation Report AS-503,
  // Apollo 8 Mission (1969),
  // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19690015314.pdf.
  // We use the reconstructed or actual values.

  // Table 2-2. Significant Event Times Summary.
  Instant const range_zero;
  Instant const s_ivb_1st_90_percent_thrust = range_zero +    530.53 * Second;
  Instant const s_ivb_1st_eco               = range_zero +    684.98 * Second;
  // Initiate S-IVB Restart Sequence and Start of Time Base 6 (T6).
  Instant const t6                          = range_zero +   9659.54 * Second;
  Instant const s_ivb_2nd_90_percent_thrust = range_zero + 10'240.02 * Second;
  Instant const s_ivb_2nd_eco               = range_zero + 10'555.51 * Second;

  // From Table 7-2. S-IVB Steady State Performance - First Burn.
  Force thrust_1st                  = 901'557 * Newton;
  Speed specific_impulse_1st        = 4'204.1 * Newton * Second / Kilogram;
  Variation<Mass> lox_flowrate_1st  = 178.16 * Kilogram / Second;
  Variation<Mass> fuel_flowrate_1st = 36.30 * Kilogram / Second;

  // From Table 7-7. S-IVB Steady State Performance - Second Burn.
  Force thrust_2nd                  = 897'548 * Newton;
  Speed specific_impulse_2nd        = 4199.2 * Newton * Second / Kilogram;
  Variation<Mass> lox_flowrate_2nd  = 177.70 * Kilogram / Second;
  Variation<Mass> fuel_flowrate_2nd = 36.01 * Kilogram / Second;

  // Table 21-5. Total Vehicle Mass, S-IVB First Burn Phase, Kilograms.
  Mass total_vehicle_at_s_ivb_1st_90_percent_thrust = 161143 * Kilogram;
  Mass total_vehicle_at_s_ivb_1st_eco               = 128095 * Kilogram;

  // Table 21-7. Total Vehicle Mass, S-IVB Second Burn Phase, Kilograms.
  Mass total_vehicle_at_s_ivb_2nd_90_percent_thrust = 126780 * Kilogram;
  Mass total_vehicle_at_s_ivb_2nd_eco               =  59285 * Kilogram;

  // An arbitrary direction, we're not testing this.
  Vector<double, Frenet<Rendering>> e_y({0, 1, 0});

  Manœuvre<World, Rendering> first_burn(
      thrust_1st,
      total_vehicle_at_s_ivb_1st_90_percent_thrust,
      specific_impulse_1st,
      e_y,
      MakeMockDynamicFrame());
  EXPECT_THAT(RelativeError(lox_flowrate_1st + fuel_flowrate_1st,
                            first_burn.mass_flow()),
              Lt(1E-4));

  first_burn.set_duration(s_ivb_1st_eco - s_ivb_1st_90_percent_thrust);
  EXPECT_THAT(
      RelativeError(total_vehicle_at_s_ivb_1st_eco, first_burn.final_mass()),
      Lt(1E-3));

  first_burn.set_initial_time(s_ivb_1st_90_percent_thrust);
  EXPECT_EQ(s_ivb_1st_eco, first_burn.final_time());

  // Accelerations from Figure 4-4. Ascent Trajectory Acceleration Comparison.
  // Final acceleration from Table 4-2. Comparison of Significant Trajectory
  // Events.
  discrete_trajectory_.Append(first_burn.initial_time(), dof_);
  EXPECT_CALL(*mock_dynamic_frame_,
              ToThisFrameAtTime(first_burn.initial_time()))
      .WillOnce(Return(rigid_motion_));
  EXPECT_CALL(*mock_dynamic_frame_,
              FrenetFrame(first_burn.initial_time(), rendering_dof_))
      .WillOnce(
          Return(Rotation<Frenet<Rendering>, Rendering>::Identity()));
  auto const first_acceleration = first_burn.acceleration(discrete_trajectory_);
  EXPECT_THAT(
      first_acceleration(first_burn.initial_time()).Norm(),
      AllOf(Gt(5 * Metre / Pow<2>(Second)), Lt(6.25 * Metre / Pow<2>(Second))));
  EXPECT_THAT(first_acceleration(range_zero + 600 * Second).Norm(),
              AllOf(Gt(6.15 * Metre / Pow<2>(Second)),
                    Lt(6.35 * Metre / Pow<2>(Second))));
  EXPECT_THAT(first_acceleration(first_burn.final_time()).Norm(),
              AllOf(Gt(7.03 * Metre / Pow<2>(Second)),
                    Lt(7.05 * Metre / Pow<2>(Second))));

  Manœuvre<World, Rendering> second_burn(
      thrust_2nd,
      total_vehicle_at_s_ivb_2nd_90_percent_thrust,
      specific_impulse_2nd,
      e_y,
      MakeMockDynamicFrame());
  EXPECT_THAT(RelativeError(lox_flowrate_2nd + fuel_flowrate_2nd,
                            second_burn.mass_flow()),
              Lt(2E-4));

  second_burn.set_duration(s_ivb_2nd_eco - s_ivb_2nd_90_percent_thrust);
  EXPECT_THAT(
      RelativeError(total_vehicle_at_s_ivb_2nd_eco, second_burn.final_mass()),
      Lt(2E-3));

  second_burn.set_initial_time(s_ivb_2nd_90_percent_thrust);
  EXPECT_EQ(s_ivb_2nd_eco, second_burn.final_time());

  // Accelerations from Figure 4-9. Injection Phase Acceleration Comparison.
  // Final acceleration from Table 4-2. Comparison of Significant Trajectory
  // Events.
  discrete_trajectory_.Append(second_burn.initial_time(), dof_);
  EXPECT_CALL(*mock_dynamic_frame_,
              ToThisFrameAtTime(second_burn.initial_time()))
      .WillOnce(Return(rigid_motion_));
  EXPECT_CALL(*mock_dynamic_frame_,
              FrenetFrame(second_burn.initial_time(), rendering_dof_))
      .WillOnce(
          Return(Rotation<Frenet<Rendering>, Rendering>::Identity()));
  auto const second_acceleration =
      second_burn.acceleration(discrete_trajectory_);
  EXPECT_THAT(second_acceleration(second_burn.initial_time()).Norm(),
              AllOf(Gt(7 * Metre / Pow<2>(Second)),
                    Lt(7.5 * Metre / Pow<2>(Second))));
  EXPECT_THAT(second_acceleration(t6 + 650 * Second).Norm(),
              AllOf(Gt(8 * Metre / Pow<2>(Second)),
                    Lt(8.02 * Metre / Pow<2>(Second))));
  EXPECT_THAT(second_acceleration(t6 + 700 * Second).Norm(),
              AllOf(Gt(8.8 * Metre / Pow<2>(Second)),
                    Lt(9 * Metre / Pow<2>(Second))));
  EXPECT_THAT(second_acceleration(t6 + 750 * Second).Norm(),
              AllOf(Gt(9.9 * Metre / Pow<2>(Second)),
                    Lt(10 * Metre / Pow<2>(Second))));
  EXPECT_THAT(second_acceleration(t6 + 850 * Second).Norm(),
              AllOf(Gt(12.97 * Metre / Pow<2>(Second)),
                    Lt(13 * Metre / Pow<2>(Second))));
  EXPECT_THAT(second_acceleration(second_burn.final_time()).Norm(),
              AllOf(Gt(15.12 * Metre / Pow<2>(Second)),
                    Lt(15.17 * Metre / Pow<2>(Second))));

  EXPECT_THAT(second_burn.Δv(),
              AllOf(Gt(3 * Kilo(Metre) / Second),
                    Lt(3.25 * Kilo(Metre) / Second)));

  // From the Apollo 8 flight journal.
  EXPECT_THAT(AbsoluteError(10'519.6 * Foot / Second, second_burn.Δv()),
              Lt(20 * Metre / Second));
}

}  // namespace ksp_plugin
}  // namespace principia
