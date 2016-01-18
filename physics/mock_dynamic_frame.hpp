
#pragma once

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock-generated-function-mockers.h"
#include "gmock/gmock-matchers.h"
#include "gmock/gmock.h"  // IWYU pragma: keep
#include "ksp_plugin/manœuvre.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

namespace geometry {
template <typename FrameTag, FrameTag frame_tag, bool frame_is_inertial> class Frame;
}  // namespace geometry
namespace physics {
template <typename FromFrame, typename ToFrame> class RigidMotion;
}  // namespace physics
namespace serialization {
class DynamicFrame;
}  // namespace serialization

using geometry::Instant;
using geometry::Vector;
using quantities::Acceleration;

namespace physics {

template<typename InertialFrame, typename ThisFrame>
class MockDynamicFrame : public DynamicFrame<InertialFrame, ThisFrame> {
 public:
  MockDynamicFrame() {}

  MOCK_CONST_METHOD1_T(ToThisFrameAtTime,
                       RigidMotion<InertialFrame, ThisFrame>(Instant const& t));
  MOCK_CONST_METHOD1_T(FromThisFrameAtTime,
                       RigidMotion<ThisFrame, InertialFrame>(Instant const& t));
  MOCK_CONST_METHOD2_T(
      GeometricAcceleration,
      Vector<Acceleration, ThisFrame>(
          Instant const& t,
          DegreesOfFreedom<ThisFrame> const& degrees_of_freedom));

  using Rot = Rotation<Frenet<ThisFrame>, ThisFrame>;

  MOCK_CONST_METHOD2_T(
      FrenetFrame,
      Rot(Instant const& t,
          DegreesOfFreedom<ThisFrame> const& degrees_of_freedom));

  MOCK_CONST_METHOD1_T(
      WriteToMessage,
      void(not_null<serialization::DynamicFrame*> const message));
};

}  // namespace physics
}  // namespace principia
