﻿
#pragma once

#include "physics/barycentric_rotating_dynamic_frame.hpp"

#include "geometry/barycentre_calculator.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_barycentric_rotating_dynamic_frame {

using geometry::Barycentre;
using geometry::Bivector;
using geometry::Displacement;
using geometry::R3x3Matrix;
using geometry::Velocity;
using geometry::Wedge;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Pow;
using quantities::Product;
using quantities::Speed;
using quantities::si::Radian;

template<typename InertialFrame, typename ThisFrame>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
BarycentricRotatingDynamicFrame(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    not_null<MassiveBody const*> const primary,
    not_null<MassiveBody const*> const secondary)
    : ephemeris_(ephemeris),
      primary_(primary),
      secondary_(secondary),
      primary_trajectory_(ephemeris_->trajectory(primary_)),
      secondary_trajectory_(ephemeris_->trajectory(secondary_)) {}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const primary_degrees_of_freedom =
      primary_trajectory_->EvaluateDegreesOfFreedom(t, &primary_hint_);
  DegreesOfFreedom<InertialFrame> const secondary_degrees_of_freedom =
      secondary_trajectory_->EvaluateDegreesOfFreedom(t, &secondary_hint_);
  DegreesOfFreedom<InertialFrame> const barycentre_degrees_of_freedom =
      Barycentre<DegreesOfFreedom<InertialFrame>, GravitationalParameter>(
          {primary_degrees_of_freedom,
           secondary_degrees_of_freedom},
          {primary_->gravitational_parameter(),
           secondary_->gravitational_parameter()});

  Rotation<InertialFrame, ThisFrame> rotation =
          Rotation<InertialFrame, ThisFrame>::Identity();
  AngularVelocity<InertialFrame> angular_velocity;
  ComputeAngularDegreesOfFreedom(primary_degrees_of_freedom,
                                 secondary_degrees_of_freedom,
                                 &rotation,
                                 &angular_velocity);

  RigidTransformation<InertialFrame, ThisFrame> const
      rigid_transformation(barycentre_degrees_of_freedom.position(),
                           ThisFrame::origin,
                           rotation.Forget());
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             angular_velocity,
             barycentre_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
void BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
WriteToMessage(not_null<serialization::DynamicFrame*> const message) const {
  auto* const extension =
      message->MutableExtension(
          serialization::BarycentricRotatingDynamicFrame::
              barycentric_rotating_dynamic_frame);
  extension->set_primary(ephemeris_->serialization_index_for_body(primary_));
  extension->set_secondary(
      ephemeris_->serialization_index_for_body(secondary_));
}

template<typename InertialFrame, typename ThisFrame>
not_null<std::unique_ptr<
    BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>>>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::BarycentricRotatingDynamicFrame const & message) {
  return std::make_unique<BarycentricRotatingDynamicFrame>(
      ephemeris,
      ephemeris->body_for_serialization_index(message.primary()),
      ephemeris->body_for_serialization_index(message.secondary()));
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, InertialFrame>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
    GravitationalAcceleration(Instant const& t,
                              Position<InertialFrame> const& q) const {
  return ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(q, t);
}

template<typename InertialFrame, typename ThisFrame>
SecondOrderRigidMotion<InertialFrame, ThisFrame>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::Motion(
    Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const primary_degrees_of_freedom =
      primary_trajectory_->EvaluateDegreesOfFreedom(t, &primary_hint_);
  DegreesOfFreedom<InertialFrame> const secondary_degrees_of_freedom =
      secondary_trajectory_->EvaluateDegreesOfFreedom(t, &secondary_hint_);

  Vector<Acceleration, InertialFrame> const primary_acceleration =
      ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(primary_, t);
  Vector<Acceleration, InertialFrame> const secondary_acceleration =
      ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(secondary_, t);

  auto const to_this_frame = ToThisFrameAtTime(t);

  // TODO(egg): TeX and reference.
  RelativeDegreesOfFreedom<InertialFrame> const primary_secondary =
      primary_degrees_of_freedom - secondary_degrees_of_freedom;
  Variation<AngularVelocity<InertialFrame>> const
      angular_acceleration_of_to_frame =
          (Wedge(primary_secondary.displacement(),
                 (primary_acceleration - secondary_acceleration)) * Radian -
           2 * to_this_frame.angular_velocity_of_to_frame() *
               InnerProduct(primary_secondary.displacement(),
                            primary_secondary.velocity())) /
          InnerProduct(primary_secondary.displacement(),
                       primary_secondary.displacement());

  Vector<Acceleration, InertialFrame> const acceleration_of_to_frame_origin =
      Barycentre<Vector<Acceleration, InertialFrame>, GravitationalParameter>(
          {primary_acceleration, secondary_acceleration},
          {primary_->gravitational_parameter(),
           secondary_->gravitational_parameter()});
  return SecondOrderRigidMotion<InertialFrame, ThisFrame>(
             to_this_frame,
             angular_acceleration_of_to_frame,
             acceleration_of_to_frame_origin);
}

template<typename InertialFrame, typename ThisFrame>
void BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
ComputeAngularDegreesOfFreedom(
    DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
    DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
    not_null<Rotation<InertialFrame, ThisFrame>*> const rotation,
    not_null<AngularVelocity<InertialFrame>*> const angular_velocity) {
  RelativeDegreesOfFreedom<InertialFrame> const reference =
      primary_degrees_of_freedom - secondary_degrees_of_freedom;
  Displacement<InertialFrame> const& reference_direction =
      reference.displacement();
  Velocity<InertialFrame> reference_normal = reference.velocity();
  reference_direction.template Orthogonalize<Speed>(&reference_normal);
  Bivector<Product<Length, Speed>, InertialFrame> const reference_binormal =
      Wedge(reference_direction, reference_normal);
  *rotation = Rotation<InertialFrame, ThisFrame>(
                  R3x3Matrix(Normalize(reference_direction).coordinates(),
                             Normalize(reference_normal).coordinates(),
                             Normalize(reference_binormal).coordinates()));
  *angular_velocity = reference_binormal * Radian /
                      InnerProduct(reference_direction, reference_direction);
}

}  // namespace internal_barycentric_rotating_dynamic_frame
}  // namespace physics
}  // namespace principia
