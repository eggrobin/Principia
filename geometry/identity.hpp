﻿
#pragma once

#include "base/mappable.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/linear_map.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/sign.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
class OrthogonalMap;

// The identity map.
template<typename FromFrame, typename ToFrame>
class Identity : public LinearMap<FromFrame, ToFrame> {
 public:
  Identity();
  ~Identity() override = default;

  Sign Determinant() const override;

  Identity<ToFrame, FromFrame> Inverse() const;

  template<typename Scalar>
  Vector<Scalar, ToFrame> operator()(
      Vector<Scalar, FromFrame> const& vector) const;

  template<typename Scalar>
  Bivector<Scalar, ToFrame> operator()(
      Bivector<Scalar, FromFrame> const& bivector) const;

  template<typename Scalar>
  Trivector<Scalar, ToFrame> operator()(
      Trivector<Scalar, FromFrame> const& trivector) const;

  template<typename T>
  typename base::Mappable<Identity, T>::type operator()(T const& t) const;

  OrthogonalMap<FromFrame, ToFrame> Forget() const;

  void WriteToMessage(not_null<serialization::LinearMap*> const message) const;
  static Identity ReadFromMessage(serialization::LinearMap const& message);

  void WriteToMessage(not_null<serialization::Identity*> const message) const;
  static Identity ReadFromMessage(serialization::Identity const& message);

 private:
  template<typename Scalar>
  R3Element<Scalar> operator()(R3Element<Scalar> const& r3_element) const;
};

template<typename FromFrame, typename ThroughFrame, typename ToFrame>
Identity<FromFrame, ToFrame> operator*(
    Identity<ThroughFrame, ToFrame> const& left,
    Identity<FromFrame, ThroughFrame> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/identity_body.hpp"
