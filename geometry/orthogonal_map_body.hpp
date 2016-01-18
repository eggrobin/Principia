﻿
#pragma once

#include "base/mappable.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "google/protobuf/extension_set.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/rotating_body.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

template <typename FromFrame, typename ToFrame> class OrthogonalMap;

template<typename FromFrame, typename ToFrame>
Sign OrthogonalMap<FromFrame, ToFrame>::Determinant() const {
  return determinant_;
}

template<typename FromFrame, typename ToFrame>
Rotation<FromFrame, ToFrame> const&
OrthogonalMap<FromFrame, ToFrame>::rotation() const {
  return rotation_;
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<ToFrame, FromFrame>
OrthogonalMap<FromFrame, ToFrame>::Inverse() const {
  return OrthogonalMap<ToFrame, FromFrame>(determinant_, rotation_.Inverse());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>(determinant_ * rotation_(vector));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>(rotation_(bivector));
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> OrthogonalMap<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return determinant_ * trivector;
}

template<typename FromFrame, typename ToFrame>
template<typename T>
typename base::Mappable<OrthogonalMap<FromFrame, ToFrame>, T>::type
OrthogonalMap<FromFrame, ToFrame>::operator()(T const& t) const {
  return base::Mappable<OrthogonalMap, T>::Do(*this, t);
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame>
OrthogonalMap<FromFrame, ToFrame>::Identity() {
  return OrthogonalMap(Sign(1), Rotation<FromFrame, ToFrame>::Identity());
}

template<typename FromFrame, typename ToFrame>
void OrthogonalMap<FromFrame, ToFrame>::WriteToMessage(
      not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(
      message->MutableExtension(serialization::OrthogonalMap::extension));
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame>
OrthogonalMap<FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::OrthogonalMap::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::OrthogonalMap::extension));
}

template<typename FromFrame, typename ToFrame>
void OrthogonalMap<FromFrame, ToFrame>::WriteToMessage(
      not_null<serialization::OrthogonalMap*> const message) const {
  determinant_.WriteToMessage(message->mutable_determinant());
  rotation_.WriteToMessage(message->mutable_rotation());
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame>
OrthogonalMap<FromFrame, ToFrame>::ReadFromMessage(
    serialization::OrthogonalMap const& message) {
  return OrthogonalMap(Sign::ReadFromMessage(message.determinant()),
                       Rotation<FromFrame, ToFrame>::ReadFromMessage(
                           message.rotation()));
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame>::OrthogonalMap(
    Sign const& determinant,
    Rotation<FromFrame, ToFrame> const& rotation)
    : determinant_(determinant),
      rotation_(rotation) {}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> operator*(
    OrthogonalMap<ThroughFrame, ToFrame> const& left,
    OrthogonalMap<FromFrame, ThroughFrame> const& right) {
  return OrthogonalMap<FromFrame, ToFrame>(
             left.determinant_ * right.determinant_,
             left.rotation_ * right.rotation_);
}

}  // namespace geometry
}  // namespace principia
