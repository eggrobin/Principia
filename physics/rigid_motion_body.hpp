﻿
#pragma once

#include "physics/rigid_motion.hpp"

#include "geometry/linear_map.hpp"

namespace principia {
namespace physics {
namespace internal_rigid_motion {

using geometry::LinearMap;

template<typename FromFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame>::RigidMotion(
    RigidTransformation<FromFrame, ToFrame> const& rigid_transformation,
    AngularVelocity<FromFrame> const& angular_velocity_of_to_frame,
    Velocity<FromFrame> const& velocity_of_to_frame_origin)
    : rigid_transformation_(rigid_transformation),
      angular_velocity_of_to_frame_(angular_velocity_of_to_frame),
      velocity_of_to_frame_origin_(velocity_of_to_frame_origin) {}

template<typename FromFrame, typename ToFrame>
RigidTransformation<FromFrame, ToFrame> const&
RigidMotion<FromFrame, ToFrame>::rigid_transformation() const {
  return rigid_transformation_;
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> const&
RigidMotion<FromFrame, ToFrame>::orthogonal_map() const {
  return rigid_transformation_.linear_map();
}

template<typename FromFrame, typename ToFrame>
AngularVelocity<FromFrame> const&
RigidMotion<FromFrame, ToFrame>::angular_velocity_of_to_frame() const {
  return angular_velocity_of_to_frame_;
}

template<typename FromFrame, typename ToFrame>
Velocity<FromFrame> const&
RigidMotion<FromFrame, ToFrame>::velocity_of_to_frame_origin() const {
  return velocity_of_to_frame_origin_;
}

template<typename FromFrame, typename ToFrame>
DegreesOfFreedom<ToFrame> RigidMotion<FromFrame, ToFrame>::operator()(
    DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const {
  return {rigid_transformation_(degrees_of_freedom.position()),
          orthogonal_map()(
              degrees_of_freedom.velocity() - velocity_of_to_frame_origin_ -
              angular_velocity_of_to_frame_ *
                  (degrees_of_freedom.position() -
                   rigid_transformation_.Inverse()(ToFrame::origin)) /
                  Radian)};
}

template<typename FromFrame, typename ToFrame>
RigidMotion<ToFrame, FromFrame>
RigidMotion<FromFrame, ToFrame>::Inverse() const {
  return RigidMotion<ToFrame, FromFrame>(
      rigid_transformation_.Inverse(),
      -orthogonal_map()(angular_velocity_of_to_frame_),
      (*this)({FromFrame::origin, Velocity<FromFrame>()}).velocity());
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame> operator*(
    RigidMotion<ThroughFrame, ToFrame> const& left,
    RigidMotion<FromFrame, ThroughFrame> const& right) {
  return RigidMotion<FromFrame, ToFrame>(
      left.rigid_transformation() * right.rigid_transformation(),
      right.angular_velocity_of_to_frame_ +
          right.orthogonal_map().Inverse()(left.angular_velocity_of_to_frame_),
      right.Inverse()(left.Inverse()(
          {ToFrame::origin, Velocity<ToFrame>()})).velocity());
}

template<typename FromFrame, typename ToFrame>
SecondOrderRigidMotion<FromFrame, ToFrame>::SecondOrderRigidMotion(
    RigidMotion<FromFrame, ToFrame> const& first_order_motion,
    Variation<AngularVelocity<FromFrame>> const&
        angular_acceleration_of_to_frame,
    Vector<Acceleration, FromFrame> const& acceleration_of_to_frame_origin)
    : first_order_motion_(first_order_motion),
      angular_acceleration_of_to_frame_(angular_acceleration_of_to_frame),
      acceleration_of_to_frame_origin_(acceleration_of_to_frame_origin) {}

template<typename FromFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame>
SecondOrderRigidMotion<FromFrame, ToFrame>::first_order_motion() const {
  return first_order_motion_;
}
template<typename FromFrame, typename ToFrame>
Variation<AngularVelocity<FromFrame>>
SecondOrderRigidMotion<FromFrame, ToFrame>::angular_acceleration_of_to_frame()
    const {
  return angular_acceleration_of_to_frame_;
}
template<typename FromFrame, typename ToFrame>
Vector<Acceleration, FromFrame>
SecondOrderRigidMotion<FromFrame, ToFrame>::acceleration_of_to_frame_origin()
    const {
  return acceleration_of_to_frame_origin_;
}

}  // namespace internal_rigid_motion
}  // namespace physics
}  // namespace principia
