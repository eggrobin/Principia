﻿#pragma once

#include <functional>

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rotation.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::AffineMap;
using geometry::AngularVelocity;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Vector;
using quantities::Acceleration;
using quantities::Length;
using quantities::si::Radian;

namespace physics {

// An arbitrary rigid transformation.  Simultaneous positions between two frames
// are always related by such a transformation.
template<typename FromFrame, typename ToFrame>
using RigidTransformation =
    AffineMap<FromFrame, ToFrame, Length, OrthogonalMap>;

// The instantaneous motion of |ToFrame| with respect to |FromFrame|.
// This is the derivative of a |RigidTransformation<FromFrame, ToFrame>|.
// In order to invert, the |RigidTransformation| is needed, and we need its
// linear part anyway, so we store it (and we forward its action on positions).
template<typename FromFrame, typename ToFrame>
class RigidMotion {
 public:
  RigidMotion(
      RigidTransformation<FromFrame, ToFrame> const& rigid_transformation,
      AngularVelocity<FromFrame> const& angular_velocity_of_to_frame,
      Velocity<FromFrame> const& velocity_of_to_frame_origin);
  ~RigidMotion() = default;

  RigidTransformation<FromFrame, ToFrame> const& rigid_transformation() const;
  // Returns |rigid_transformation().linear_map()|.
  OrthogonalMap<FromFrame, ToFrame> const& orthogonal_map() const;
  AngularVelocity<FromFrame> const& angular_velocity_of_to_frame() const;
  Velocity<FromFrame> const& velocity_of_to_frame_origin() const;

  DegreesOfFreedom<ToFrame> operator()(
      DegreesOfFreedom<FromFrame> const& degrees_of_freedom) const;

  RigidMotion<ToFrame, FromFrame> Inverse() const;

 private:
  RigidTransformation<FromFrame, ToFrame> const rigid_transformation_;
  // d/dt rigid_transformation⁻¹(basis of ToFrame). The positively oriented
  // orthogonal bases of |FromFrame| are acted upon faithfully and transitively
  // by SO(FromFrame), so this lies in the tangent space, i.e., the Lie algebra
  // 𝖘𝔬(FromFrame) ≅ FromFrame ∧ FromFrame.
  AngularVelocity<FromFrame> angular_velocity_of_to_frame_;
  // d/dt rigid_transformation⁻¹(ToFrame::origin).
  Velocity<FromFrame> velocity_of_to_frame_origin_;

  template<typename From, typename Through, typename To>
  friend RigidMotion<From, To> operator*(
      RigidMotion<Through, To> const& left,
      RigidMotion<From, Through> const& right);
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
RigidMotion<FromFrame, ToFrame> operator*(
    RigidMotion<ThroughFrame, ToFrame> const& left,
    RigidMotion<FromFrame, ThroughFrame> const& right);

}  // namespace physics
}  // namespace principia

#include "physics/rigid_motion_body.hpp"
