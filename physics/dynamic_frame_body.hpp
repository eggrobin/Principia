﻿
#pragma once

#include "physics/barycentric_rotating_dynamic_frame.hpp"
#include "physics/body_centered_non_rotating_dynamic_frame.hpp"
#include "physics/dynamic_frame.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_dynamic_frame {

using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::Displacement;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::R3x3Matrix;
using geometry::Velocity;
using geometry::Wedge;
using quantities::Pow;
using quantities::Sqrt;
using quantities::si::Radian;

template<typename InertialFrame, typename ThisFrame>
RigidMotion<ThisFrame, InertialFrame>
DynamicFrame<InertialFrame, ThisFrame>::FromThisFrameAtTime(
    Instant const& t) const {
  return ToThisFrameAtTime(t).Inverse();
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
DynamicFrame<InertialFrame, ThisFrame>::GeometricAcceleration(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  auto const motion = Motion(t);
  auto const to_this_frame = motion.first_order_motion();
  auto const from_this_frame = to_this_frame.Inverse();

  // Beware, we want the angular velocity of ThisFrame as seen in the
  // InertialFrame, but pushed to ThisFrame.  Otherwise the sign is wrong.
  AngularVelocity<InertialFrame> const Ω_inertial =
      to_this_frame.angular_velocity_of_to_frame();
  AngularVelocity<ThisFrame> const Ω =
      to_this_frame.orthogonal_map()(Ω_inertial);

  Variation<AngularVelocity<InertialFrame>> const dΩ_over_dt_inertial =
      motion.angular_acceleration_of_to_frame();
  Variation<AngularVelocity<ThisFrame>> const dΩ_over_dt =
      to_this_frame.orthogonal_map()(dΩ_over_dt_inertial);

  Vector<Acceleration, ThisFrame> const gravitational_acceleration_at_point =
      to_this_frame.orthogonal_map()(
          GravitationalAcceleration(t,
                                    from_this_frame.rigid_transformation()(
                                        degrees_of_freedom.position())));

  Vector<Acceleration, ThisFrame> const linear_acceleration =
      -to_this_frame.orthogonal_map()(motion.acceleration_of_to_frame_origin());

  Displacement<ThisFrame> const r =
      degrees_of_freedom.position() - ThisFrame::origin;
  Vector<Acceleration, ThisFrame> const coriolis_acceleration_at_point =
      -2 * Ω * degrees_of_freedom.velocity() / Radian;
  Vector<Acceleration, ThisFrame> const centrifugal_acceleration_at_point =
      -Ω * (Ω * r) / Pow<2>(Radian);
  Vector<Acceleration, ThisFrame> const euler_acceleration_at_point =
      -dΩ_over_dt * r / Radian;

  Vector<Acceleration, ThisFrame> const fictitious_acceleration =
      linear_acceleration +
      coriolis_acceleration_at_point +
      centrifugal_acceleration_at_point +
      euler_acceleration_at_point;
  return gravitational_acceleration_at_point + fictitious_acceleration;
}

template<typename InertialFrame, typename ThisFrame>
Rotation<Frenet<ThisFrame>, ThisFrame>
DynamicFrame<InertialFrame, ThisFrame>::FrenetFrame(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  Velocity<ThisFrame> const& velocity = degrees_of_freedom.velocity();
  Vector<Acceleration, ThisFrame> const acceleration =
      GeometricAcceleration(t, degrees_of_freedom);
  Vector<Acceleration, ThisFrame> normal_acceleration = acceleration;
  velocity.template Orthogonalize<Acceleration>(&normal_acceleration);
  Vector<double, ThisFrame> tangent = Normalize(velocity);
  Vector<double, ThisFrame> normal = Normalize(normal_acceleration);
  Bivector<double, ThisFrame> binormal = Wedge(tangent, normal);
  // Maps |tangent| to {1, 0, 0}, |normal| to {0, 1, 0}, and |binormal| to
  // {0, 0, 1}.
  return Rotation<Frenet<ThisFrame>, ThisFrame>(
      R3x3Matrix(tangent.coordinates(),
                 normal.coordinates(),
                 binormal.coordinates()).Transpose());
}

template<typename InertialFrame, typename ThisFrame>
std::unique_ptr<DynamicFrame<InertialFrame, ThisFrame>>
DynamicFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::DynamicFrame const& message) {
  std::unique_ptr<DynamicFrame> result;
  int extensions_found = 0;
  // NOTE(egg): the |static_cast|ing below is needed on MSVC, because the silly
  // compiler doesn't see the |operator std::unique_ptr<DynamicFrame>() &&|.
  if (message.HasExtension(serialization::BarycentricRotatingDynamicFrame::
                               barycentric_rotating_dynamic_frame)) {
    ++extensions_found;
    result = static_cast<not_null<std::unique_ptr<DynamicFrame>>>(
        BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
            ReadFromMessage(ephemeris,
                            message.GetExtension(
                                serialization::BarycentricRotatingDynamicFrame::
                                    barycentric_rotating_dynamic_frame)));
  }
  if (message.HasExtension(serialization::BodyCentredNonRotatingDynamicFrame::
                               body_centred_non_rotating_dynamic_frame)) {
    ++extensions_found;
    result = static_cast<not_null<std::unique_ptr<DynamicFrame>>>(
        BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
            ReadFromMessage(
                ephemeris,
                message.GetExtension(
                    serialization::BodyCentredNonRotatingDynamicFrame::
                        body_centred_non_rotating_dynamic_frame)));
  }
  CHECK_LE(extensions_found, 1) << message.DebugString();
  // For pre-Brouwer compatibility, return a null pointer if no extension is
  // found.
  return result;
}

}  // namespace internal_dynamic_frame
}  // namespace physics
}  // namespace principia
