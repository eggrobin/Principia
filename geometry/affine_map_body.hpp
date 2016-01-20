﻿
#pragma once

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/rotation.hpp"
#include "mathematica/mathematica.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

// The map is represented as x ↦ linear_map(x - from_origin) + to_origin.  This
// numerically better behaved than x ↦ linear_map(x) + translation with
// translation = to_origin - linear_map(from_origin).
template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::AffineMap(
    Point<FromVector> const& from_origin,
    Point<ToVector> const& to_origin,
    LinearMap<FromFrame, ToFrame> const& linear_map)
    : from_origin_(from_origin),
      to_origin_(to_origin),
      linear_map_(linear_map) {}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<ToFrame, FromFrame, Scalar, LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::Inverse() const {
  return AffineMap<ToFrame, FromFrame, Scalar, LinearMap>(
      to_origin_ /*from_origin*/,
      from_origin_ /*to_origin*/,
      linear_map_.Inverse() /*linear_map*/);
}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
Point<typename AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::ToVector>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::operator()(
    Point<FromVector> const& point) const {
  return Point<
      typename AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::ToVector>(
          linear_map_(point - from_origin_) + to_origin_);
}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::Identity() {
  return AffineMap();
}

template <typename FromFrame,
          typename ToFrame,
          typename Scalar,
          template <typename, typename> class LinearMap>
LinearMap<FromFrame, ToFrame> const&
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::linear_map() const {
  return linear_map_;
}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
void AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::WriteToMessage(
    not_null<serialization::AffineMap*> const message) const {
  FromFrame::WriteToMessage(message->mutable_from_frame());
  ToFrame::WriteToMessage(message->mutable_to_frame());
  from_origin_.WriteToMessage(message->mutable_from_origin());
  to_origin_.WriteToMessage(message->mutable_to_origin());
  linear_map_.WriteToMessage(message->mutable_linear_map());
}

template<typename FromFrame, typename ToFrame, typename Scalar,
         template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap>::ReadFromMessage(
    serialization::AffineMap const& message) {
  FromFrame::ReadFromMessage(message.from_frame());
  ToFrame::ReadFromMessage(message.to_frame());
  return AffineMap(Point<FromVector>::ReadFromMessage(message.from_origin()),
                   Point<ToVector>::ReadFromMessage(message.to_origin()),
                   LinearMap<FromFrame, ToFrame>::ReadFromMessage(
                       message.linear_map()));
}

template<typename FromFrame, typename ThroughFrame, typename ToFrame,
         typename Scalar, template<typename, typename> class LinearMap>
AffineMap<FromFrame, ToFrame, Scalar, LinearMap> operator*(
    AffineMap<ThroughFrame, ToFrame, Scalar, LinearMap> const& left,
    AffineMap<FromFrame, ThroughFrame, Scalar, LinearMap> const& right) {
  return AffineMap<FromFrame, ToFrame, Scalar, LinearMap>(
      right.from_origin_ /*from_origin*/,
      left.to_origin_ +
          left.linear_map_(right.to_origin_ - left.from_origin_) /*to_origin*/,
      left.linear_map_ * right.linear_map_ /*linear_map*/);
}

}  // namespace geometry
}  // namespace principia
