﻿
#pragma once

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::LinearMap;

namespace physics {

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

}  // namespace physics
}  // namespace principia
