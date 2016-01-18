﻿#pragma once

// We use ostream for logging purposes.
#include <iostream>  // NOLINT(readability/streams)
#include <string>

#include "base/not_null.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

namespace serialization {
class R3Element;
}  // namespace serialization

using base::not_null;
using quantities::Angle;

namespace geometry {

template<typename Scalar>
struct SphericalCoordinates;

// An |R3Element<Scalar>| is an element of Scalar³. |Scalar| should be a vector
// space over ℝ, represented by |double|. |R3Element| is the underlying data
// type for more advanced strongly typed structures suchas |Multivector|.
template<typename Scalar>
struct R3Element {
 public:
  R3Element();
  R3Element(Scalar const& x, Scalar const& y, Scalar const& z);

  Scalar&       operator[](int const index);
  Scalar const& operator[](int const index) const;

  R3Element& operator+=(R3Element const& right);
  R3Element& operator-=(R3Element const& right);
  R3Element& operator*=(double const right);
  R3Element& operator/=(double const right);

  Scalar Norm() const;

  // Uses the x-y plane as the equator, the x-axis as the reference direction on
  // the equator, and the z-axis as the north pole.
  SphericalCoordinates<Scalar> ToSpherical() const;

  // Modifies |*r3_element| so as to make it orthogonal to |*this|, using the
  // modified Gram-Schmidt algorithm.  Fails if |*this| is zero.
  template<typename S>
  void Orthogonalize(not_null<R3Element<S>*> const r3_element) const;

  void WriteToMessage(not_null<serialization::R3Element*> const message) const;
  static R3Element ReadFromMessage(serialization::R3Element const& message);

  Scalar x;
  Scalar y;
  Scalar z;
};

template<typename Scalar>
struct SphericalCoordinates {
  // Default, but prevents aggregate initialization of |SphericalCoordinates| to
  // obviate confusion over the order of |latitude| and |longitude|.
  SphericalCoordinates();

  // Uses the x-y plane as the equator, the x-axis as the reference direction on
  // the equator, and the z-axis as the north pole.
  R3Element<Scalar> ToCartesian();

  Scalar radius;
  Angle latitude;
  Angle longitude;
};

template<typename Scalar>
SphericalCoordinates<Scalar> RadiusLatitudeLongitude(Scalar const& radius,
                                                     Angle const& latitude,
                                                     Angle const& longitude);

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& right);

template<typename Scalar>
R3Element<Scalar> operator+(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator-(R3Element<Scalar> const& left,
                            R3Element<Scalar> const& right);

template<typename Scalar>
R3Element<Scalar> operator*(double const left,
                            R3Element<Scalar> const& right);
template<typename Scalar>
R3Element<Scalar> operator*(R3Element<Scalar> const& left,
                            double const right);
template<typename Scalar>
R3Element<Scalar> operator/(R3Element<Scalar> const& left,
                            double const right);

// Dimensionful multiplication |LScalar * R3Element<RScalar>| is the tensor
// product LScalar ⊗ Scalar³. Since LScalar ⊗ Scalar³ ≅ (LScalar ⊗ Scalar)³,
// the result is an R3Element<Product<LScalar, RScalar>>.
// The special case where one of the scalars is |double| is handled separately
// above in order to allow implicit conversions to |double|.
template<typename LDimension, typename RScalar>
R3Element<quantities::Product<quantities::Quantity<LDimension>, RScalar>>
operator*(quantities::Quantity<LDimension> const& left,
          R3Element<RScalar> const& right);
template<typename LScalar, typename RDimension>
R3Element<quantities::Product<LScalar, quantities::Quantity<RDimension>>>
operator*(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right);
template<typename LScalar, typename RDimension>
R3Element<quantities::Quotient<LScalar, quantities::Quantity<RDimension>>>
operator/(R3Element<LScalar> const& left,
          quantities::Quantity<RDimension> const& right);

template<typename Scalar>
bool operator==(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right);
template<typename Scalar>
bool operator!=(R3Element<Scalar> const& left,
                R3Element<Scalar> const& right);

template<typename Scalar>
R3Element<double> Normalize(R3Element<Scalar> const& r3_element);

template<typename Scalar>
std::string DebugString(R3Element<Scalar> const& r3_element);

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         R3Element<Scalar> const& r3_element);

template<typename LScalar, typename RScalar>
R3Element<quantities::Product<LScalar, RScalar>> Cross(
    R3Element<LScalar> const& left,
    R3Element<RScalar> const& right);

template<typename LScalar, typename RScalar>
quantities::Product<LScalar, RScalar> Dot(R3Element<LScalar> const& left,
                                          R3Element<RScalar> const& right);

}  // namespace geometry
}  // namespace principia

#include "geometry/r3_element_body.hpp"
