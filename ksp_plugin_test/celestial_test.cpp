﻿
#include <memory>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest-death-test.h"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Second;

namespace ksp_plugin {

class CelestialTest : public testing::Test {
 protected:
  CelestialTest()
      : body_(1 * Kilogram),
        celestial_(make_not_null_unique<Celestial>(&body_)),
        trajectory_(1 * Second, 1 * Metre) {}

  MassiveBody body_;
  not_null<std::unique_ptr<Celestial>> celestial_;
  ContinuousTrajectory<Barycentric> trajectory_;
};

using CelestialDeathTest = CelestialTest;

TEST_F(CelestialDeathTest, Uninitialized) {
  EXPECT_DEATH({celestial_->trajectory();}, "is_initialized");
  EXPECT_DEATH({celestial_->current_time_hint();}, "is_initialized");
  EXPECT_DEATH({celestial_->current_degrees_of_freedom(Instant());},
               "is_initialized");
  EXPECT_DEATH({celestial_->current_position(Instant());}, "is_initialized");
  EXPECT_DEATH({celestial_->current_velocity(Instant());}, "is_initialized");
}

TEST_F(CelestialDeathTest, OverlyInitialized) {
  EXPECT_DEATH({
    celestial_->set_trajectory(&trajectory_);
    celestial_->set_trajectory(&trajectory_);
  },
  "!is_initialized");
}

TEST_F(CelestialTest, Initialization) {
  EXPECT_FALSE(celestial_->is_initialized());
  ContinuousTrajectory<Barycentric> trajectory(1 * Second, 1 * Metre);
  celestial_->set_trajectory(&trajectory);
  EXPECT_TRUE(celestial_->is_initialized());
}

}  // namespace ksp_plugin
}  // namespace principia
