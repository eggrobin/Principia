﻿
#ifndef PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
#define PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_

#include <memory>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "physics/body.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/forkable.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"
#include "tools/generate_configuration.hpp"

namespace principia {

namespace physics {
template <typename Frame> class DegreesOfFreedom;
template <typename Frame> class Ephemeris;
}  // namespace physics
namespace serialization {
class DynamicFrame;
}  // namespace serialization

using geometry::Rotation;

namespace physics {

// The Frenet frame of a free fall trajectory in |Frame|.
// TODO(egg): this should actually depend on its template parameter somehow.
template<typename Frame>
using Frenet = geometry::Frame<serialization::Frame::PhysicsTag,
                               serialization::Frame::FRENET,
                               false /*frame_is_inertial*/>;

// The definition of a reference frame |ThisFrame| in arbitrary motion with
// respect to the inertial reference frame |InertialFrame|.
template<typename InertialFrame, typename ThisFrame>
class DynamicFrame {
  static_assert(InertialFrame::is_inertial, "InertialFrame must be inertial");

 public:
  virtual ~DynamicFrame() = default;
  virtual RigidMotion<InertialFrame, ThisFrame> ToThisFrameAtTime(
      Instant const& t) const = 0;
  virtual RigidMotion<ThisFrame, InertialFrame> FromThisFrameAtTime(
      Instant const& t) const = 0;

  // The acceleration due to the non-inertial motion of |ThisFrame| and gravity.
  // A particle in free fall follows a trajectory whose second derivative
  // is |GeometricAcceleration|.
  virtual Vector<Acceleration, ThisFrame> GeometricAcceleration(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const = 0;

  // The definition of the Frenet frame of a free fall trajectory in |ThisFrame|
  // with the given |degrees_of_freedom| at instant |t|.
  virtual Rotation<Frenet<ThisFrame>, ThisFrame> FrenetFrame(
      Instant const& t,
      DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const;

  virtual void WriteToMessage(
      not_null<serialization::DynamicFrame*> const message) const = 0;

  // Dispatches to one of the subclasses depending on the contents of the
  // message.  Returns |nullptr| if no dynamic frame extension is found.
  static std::unique_ptr<DynamicFrame>
      ReadFromMessage(not_null<Ephemeris<InertialFrame> const*> const ephemeris,
                      serialization::DynamicFrame const& message);
};

}  // namespace physics
}  // namespace principia

#include "physics/dynamic_frame_body.hpp"

#endif  // PRINCIPIA_PHYSICS_DYNAMIC_FRAME_HPP_
