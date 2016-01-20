﻿
#pragma once

#include "base/mappable.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "google/protobuf/extension_set.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/root_finders.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/rotating_body.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
Identity<FromFrame, ToFrame>::Identity() {}

template<typename FromFrame, typename ToFrame>
Sign Identity<FromFrame, ToFrame>::Determinant() const {
  return Sign(1);
}

template<typename FromFrame, typename ToFrame>
Identity<ToFrame, FromFrame> Identity<FromFrame, ToFrame>::Inverse() const {
  return Identity<ToFrame, FromFrame>();
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Vector<Scalar, ToFrame> Identity<FromFrame, ToFrame>::operator()(
    Vector<Scalar, FromFrame> const& vector) const {
  return Vector<Scalar, ToFrame>(vector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Bivector<Scalar, ToFrame> Identity<FromFrame, ToFrame>::operator()(
    Bivector<Scalar, FromFrame> const& bivector) const {
  return Bivector<Scalar, ToFrame>(bivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename Scalar>
Trivector<Scalar, ToFrame> Identity<FromFrame, ToFrame>::operator()(
    Trivector<Scalar, FromFrame> const& trivector) const {
  return Trivector<Scalar, ToFrame>(trivector.coordinates());
}

template<typename FromFrame, typename ToFrame>
template<typename T>
typename base::Mappable<Identity<FromFrame, ToFrame>, T>::type
Identity<FromFrame, ToFrame>::operator()(T const& t) const {
  return base::Mappable<Identity, T>::Do(*this, t);
}

template<typename FromFrame, typename ToFrame>
OrthogonalMap<FromFrame, ToFrame> Identity<FromFrame, ToFrame>::Forget() const {
  return OrthogonalMap<FromFrame, ToFrame>(
      Determinant(),
      Rotation<FromFrame, ToFrame>::Identity());
}

template<typename FromFrame, typename ToFrame>
void Identity<FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) const {
  LinearMap<FromFrame, ToFrame>::WriteToMessage(message);
  WriteToMessage(message->MutableExtension(serialization::Identity::extension));
}

template<typename FromFrame, typename ToFrame>
Identity<FromFrame, ToFrame> Identity<FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  LinearMap<FromFrame, ToFrame>::ReadFromMessage(message);
  CHECK(message.HasExtension(serialization::Identity::extension));
  return ReadFromMessage(
      message.GetExtension(serialization::Identity::extension));
}

template<typename FromFrame, typename ToFrame>
void Identity<FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::Identity*> const message) const {}

template<typename FromFrame, typename ToFrame>
Identity<FromFrame, ToFrame> Identity<FromFrame, ToFrame>::ReadFromMessage(
    serialization::Identity const& message) {
  return Identity();
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Identity<FromFrame, ToFrame> operator*(
    Identity<ThroughFrame, ToFrame> const& left,
    Identity<FromFrame, ThroughFrame> const& right) {
  return Identity<FromFrame, ToFrame>();
}

}  // namespace geometry
}  // namespace principia
