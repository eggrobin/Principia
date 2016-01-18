﻿#pragma once

#include <iosfwd>

#include "base/not_null.hpp"
#include "geometry/r3_element.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace serialization {
class Quaternion;
}  // namespace serialization
}  // namespace principia

namespace principia {
namespace geometry {

// An element of the skew field of quaternions ℍ (where ℝ is modeled by
// |double|).
class Quaternion {
 public:
  Quaternion();
  explicit Quaternion(double const real_part);
  Quaternion(double const real_part,
             R3Element<double> const& imaginary_part);

  double real_part() const;
  R3Element<double> const& imaginary_part() const;

  Quaternion Conjugate() const;
  Quaternion Inverse() const;

  Quaternion& operator+=(Quaternion const& right);
  Quaternion& operator-=(Quaternion const& right);
  Quaternion& operator*=(Quaternion const& right);
  Quaternion& operator/=(Quaternion const& right);

  Quaternion& operator*=(double const right);
  Quaternion& operator/=(double const right);

  void WriteToMessage(not_null<serialization::Quaternion*> const message) const;
  static Quaternion ReadFromMessage(serialization::Quaternion const& message);

 private:
  double real_part_;
  R3Element<double> imaginary_part_;
};

bool operator==(Quaternion const& left, Quaternion const& right);
bool operator!=(Quaternion const& left, Quaternion const& right);

Quaternion operator+(Quaternion const& right);
Quaternion operator-(Quaternion const& right);

Quaternion operator+(Quaternion const& left, Quaternion const& right);
Quaternion operator-(Quaternion const& left, Quaternion const& right);
Quaternion operator*(Quaternion const& left, Quaternion const& right);
Quaternion operator/(Quaternion const& left, Quaternion const& right);

Quaternion operator*(double const left,
                     Quaternion const& right);
Quaternion operator*(Quaternion const& left,
                     double const right);
Quaternion operator/(Quaternion const& left,
                     double const right);

std::ostream& operator<<(std::ostream& out, Quaternion const& quaternion);

}  // namespace geometry
}  // namespace principia

#include "geometry/quaternion_body.hpp"
