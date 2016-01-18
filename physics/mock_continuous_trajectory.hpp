#pragma once

#include <vector>

#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock-generated-function-mockers.h"
#include "gmock/gmock.h"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"

namespace principia {
namespace physics {

template<typename Frame>
class MockContinuousTrajectory : public ContinuousTrajectory<Frame> {
 public:
  MockContinuousTrajectory() : ContinuousTrajectory<Frame>() {}

  MOCK_CONST_METHOD2_T(
      EvaluateDegreesOfFreedom,
      DegreesOfFreedom<Frame>(
          Instant const& time,
          typename ContinuousTrajectory<Frame>::Hint* const hint));
};

}  // namespace physics
}  // namespace principia
