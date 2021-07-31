
#pragma once

#include "physics/body_centred_line_of_sight_dynamic_frame.hpp"

#include <algorithm>
#include <utility>

#include "geometry/named_quantities.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_body_centred_line_of_sight_dynamic_frame {

using base::dynamic_cast_not_null;
using geometry::Bivector;
using geometry::Displacement;
using geometry::Frame;
using geometry::OrthogonalMap;
using geometry::R3x3Matrix;
using geometry::Velocity;
using geometry::Wedge;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Pow;
using quantities::Product;
using quantities::Speed;
using quantities::Variation;
using quantities::si::Radian;

template<typename InertialFrame, typename ThisFrame>
BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::
BodyCentredLineOfSightDynamicFrame(
    not_null<Ephemeris<InertialFrame> const*> ephemeris,
    not_null<MassiveBody const*> const primary,
    not_null<RotatingBody<InertialFrame> const*> secondary)
    : ephemeris_(std::move(ephemeris)),
      primary_(primary),
      secondary_(std::move(secondary)),
      compute_gravitational_acceleration_on_primary_(
          [this](Position<InertialFrame> const& position, Instant const& t) {
            return ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(
                primary_, t);
          }),
      primary_trajectory_(
          [&t = *ephemeris_->trajectory(primary_)]() -> auto& { return t; }),
      secondary_trajectory_(ephemeris_->trajectory(secondary_)) {}

template<typename InertialFrame, typename ThisFrame>
BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::
BodyCentredLineOfSightDynamicFrame(
    not_null<Ephemeris<InertialFrame> const*> ephemeris,
    std::function<Trajectory<InertialFrame> const&()> primary_trajectory,
    not_null<MassiveBody const*> secondary)
    : ephemeris_(std::move(ephemeris)),
      primary_(nullptr),
      secondary_(std::move(secondary)),
      compute_gravitational_acceleration_on_primary_(
          [this](Position<InertialFrame> const& position, Instant const& t) {
            return ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(
                position, t);
          }),
      primary_trajectory_(std::move(primary_trajectory)),
      secondary_trajectory_(ephemeris_->trajectory(secondary_)) {}

template<typename InertialFrame, typename ThisFrame>
not_null<MassiveBody const*>
BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::primary()
    const {
  return primary_;
}

template<typename InertialFrame, typename ThisFrame>
not_null<RotatingBody<InertialFrame> const*>
BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::secondary()
    const {
  return secondary_;
}

template<typename InertialFrame, typename ThisFrame>
Instant BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::t_min()
    const {
  return std::max(primary_trajectory_().t_min(),
                  secondary_trajectory_->t_min());
}

template<typename InertialFrame, typename ThisFrame>
Instant BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::t_max()
    const {
  return std::min(primary_trajectory_().t_max(),
                  secondary_trajectory_->t_max());
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::
    ToThisFrameAtTime(Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const primary_degrees_of_freedom =
      primary_trajectory_().EvaluateDegreesOfFreedom(t);
  DegreesOfFreedom<InertialFrame> const secondary_degrees_of_freedom =
      secondary_trajectory_->EvaluateDegreesOfFreedom(t);

  Rotation<InertialFrame, ThisFrame> rotation =
      Rotation<InertialFrame, ThisFrame>::Identity();
  AngularVelocity<InertialFrame> angular_velocity;
  ComputeAngularDegreesOfFreedom(primary_degrees_of_freedom,
                                 secondary_degrees_of_freedom,
                                 secondary_->polar_axis(),
                                 rotation,
                                 angular_velocity);

  RigidTransformation<InertialFrame, ThisFrame> const
      rigid_transformation(primary_degrees_of_freedom.position(),
                           ThisFrame::origin,
                           rotation.template Forget<OrthogonalMap>());
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             angular_velocity,
             primary_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
void BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::
WriteToMessage(not_null<serialization::DynamicFrame*> const message) const {
  auto* const extension =
      message->MutableExtension(
          serialization::BodyCentredLineOfSightDynamicFrame::extension);
  extension->set_primary(ephemeris_->serialization_index_for_body(primary_));
  extension->set_secondary(
      ephemeris_->serialization_index_for_body(secondary_));
}

template<typename InertialFrame, typename ThisFrame>
not_null<std::unique_ptr<
    BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>>>
BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::BodyCentredLineOfSightDynamicFrame const& message) {
  return std::make_unique<BodyCentredLineOfSightDynamicFrame>(
      ephemeris,
      ephemeris->body_for_serialization_index(message.primary()),
      dynamic_cast_not_null<RotatingBody<InertialFrame> const*>(
          ephemeris->body_for_serialization_index(message.secondary())));
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, InertialFrame>
BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::
GravitationalAcceleration(Instant const& t,
                          Position<InertialFrame> const& q) const {
  return ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(q, t);
}

template<typename InertialFrame, typename ThisFrame>
AcceleratedRigidMotion<InertialFrame, ThisFrame>
BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::
MotionOfThisFrame(Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const primary_degrees_of_freedom =
      primary_trajectory_().EvaluateDegreesOfFreedom(t);
  DegreesOfFreedom<InertialFrame> const secondary_degrees_of_freedom =
      secondary_trajectory_->EvaluateDegreesOfFreedom(t);

  // TODO(egg): eventually we want to add the intrinsic acceleration here.
  Vector<Acceleration, InertialFrame> const primary_acceleration =
      compute_gravitational_acceleration_on_primary_(
          primary_degrees_of_freedom.position(), t);

  Vector<Acceleration, InertialFrame> const secondary_acceleration =
      ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(secondary_, t);

  auto const to_this_frame = ToThisFrameAtTime(t);

  // TODO(egg): TeX and reference.
  RelativeDegreesOfFreedom<InertialFrame> const secondary_primary =
      secondary_degrees_of_freedom - primary_degrees_of_freedom;
  Displacement<InertialFrame> const& r = secondary_primary.displacement();
  Velocity<InertialFrame> const& ṙ = secondary_primary.velocity();
  Vector<Acceleration, InertialFrame> const r̈ =
      secondary_acceleration - primary_acceleration;
  AngularVelocity<InertialFrame> const& ω =
      to_this_frame.template angular_velocity_of<ThisFrame>();
  Variation<AngularVelocity<InertialFrame>> const
      angular_acceleration_of_to_frame =
          (Wedge(r, r̈) * Radian - 2 * ω * InnerProduct(r, ṙ)) / r.Norm²();

  Vector<Acceleration, InertialFrame> const& acceleration_of_to_frame_origin =
      primary_acceleration;
  return AcceleratedRigidMotion<InertialFrame, ThisFrame>(
             to_this_frame,
             angular_acceleration_of_to_frame,
             acceleration_of_to_frame_origin);
}

template<typename InertialFrame, typename ThisFrame>
void BodyCentredLineOfSightDynamicFrame<InertialFrame, ThisFrame>::
    ComputeAngularDegreesOfFreedom(
        DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
        DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
        Vector<double, InertialFrame> const& polar_axis,
        Rotation<InertialFrame, ThisFrame>& rotation,
        AngularVelocity<InertialFrame>& angular_velocity) {
  RelativeDegreesOfFreedom<InertialFrame> const dof =
      secondary_degrees_of_freedom - primary_degrees_of_freedom;
  Displacement<InertialFrame> const& line_of_sight = dof.displacement();
  Vector<double, InertialFrame> const local_north =
      Normalize(polar_axis.OrthogonalizationAgainst(line_of_sight));
  rotation = Rotation<InertialFrame, ThisFrame>(
      local_north,
      Normalize(Wedge(line_of_sight, local_north)),
      Normalize(line_of_sight));
  // The radial motion of primary as seen from secondary does not result in any
  // rotation of this frame.
  // The longitudinal motion (in right ascension) results in rotation about the
  // celestial pole. The latitudinal motion (in declination) results in rotation
  // about the local East direction.
  // First drop the radial component.
  Velocity<InertialFrame> const horizontal_velocity =
      dof.velocity().OrthogonalizationAgainst(line_of_sight);
  Velocity<InertialFrame> latitudinal_velocity =
      InnerProduct(horizontal_velocity, local_north) * local_north;
  Velocity<InertialFrame> longitudinal_velocity =
      horizontal_velocity - latitudinal_velocity;
  Displacement<InertialFrame> parallel_radius =
      line_of_sight.OrthogonalizationAgainst(polar_axis);
  angular_velocity = Wedge(line_of_sight, latitudinal_velocity) * Radian /
                         line_of_sight.Norm²() +
                     Wedge(parallel_radius, longitudinal_velocity) * Radian /
                         parallel_radius.Norm²();
}

}  // namespace internal_body_centred_line_of_sight_dynamic_frame
}  // namespace physics
}  // namespace principia
