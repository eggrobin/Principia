﻿
#pragma once

#include <ostream>

#include "base/not_null.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/r3_element.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

class Quaternion;

inline Quaternion::Quaternion() : real_part_(0) {}

inline Quaternion::Quaternion(double const real_part)
    : real_part_(real_part) {}

inline Quaternion::Quaternion(double const real_part,
                              R3Element<double> const& imaginary_part)
    : real_part_(real_part),
      imaginary_part_(imaginary_part) {}

inline double Quaternion::real_part() const {
  return real_part_;
}

inline R3Element<double> const&
Quaternion::imaginary_part() const {
  return imaginary_part_;
}

inline Quaternion Quaternion::Conjugate() const {
  return Quaternion(real_part_, -imaginary_part_);
}

inline Quaternion Quaternion::Inverse() const {
  return Conjugate() /
      (real_part_ * real_part_ + Dot(imaginary_part_, imaginary_part_));
}

inline Quaternion& Quaternion::operator+=(Quaternion const& right) {
  real_part_ += right.real_part_;
  imaginary_part_ += right.imaginary_part_;
  return *this;
}

inline Quaternion& Quaternion::operator-=(Quaternion const& right)  {
  real_part_ -= right.real_part_;
  imaginary_part_ -= right.imaginary_part_;
  return *this;
}

inline Quaternion& Quaternion::operator*=(Quaternion const& right) {
  // TODO(phl): Can this be optimized?
  return *this = *this * right;
}

inline Quaternion& Quaternion::operator/=(Quaternion const& right) {
  // TODO(phl): Can this be optimized?
  return *this = *this / right;
}

inline Quaternion& Quaternion::operator*=(double const right) {
  real_part_ *= right;
  imaginary_part_ *= right;
  return *this;
}

inline Quaternion& Quaternion::operator/=(double const right) {
  real_part_ /= right;
  imaginary_part_ /= right;
  return *this;
}

inline void Quaternion::WriteToMessage(
    not_null<serialization::Quaternion*> const message) const {
  message->set_real_part(real_part_);
  imaginary_part_.WriteToMessage(message->mutable_imaginary_part());
}

inline Quaternion Quaternion::ReadFromMessage(
    serialization::Quaternion const& message) {
  return Quaternion(message.real_part(),
                    R3Element<double>::ReadFromMessage(
                        message.imaginary_part()));
}

inline bool operator==(Quaternion const& left, Quaternion const& right) {
  return left.real_part() == right.real_part() &&
         left.imaginary_part() == left.imaginary_part();
}

inline bool operator!=(Quaternion const& left, Quaternion const& right) {
  return left.real_part() != right.real_part() ||
         left.imaginary_part() != left.imaginary_part();
}

inline Quaternion operator+(Quaternion const& right) {
  return right;
}

inline Quaternion operator-(Quaternion const& right) {
  return Quaternion(-right.real_part(), -right.imaginary_part());
}

inline Quaternion operator+(Quaternion const& left, Quaternion const& right) {
  return Quaternion(left.real_part() + right.real_part(),
                    left.imaginary_part() + right.imaginary_part());
}

inline Quaternion operator-(Quaternion const& left, Quaternion const& right) {
  return Quaternion(left.real_part() - right.real_part(),
                    left.imaginary_part() - right.imaginary_part());
}

inline Quaternion operator*(Quaternion const& left, Quaternion const& right) {
  return Quaternion(left.real_part() * right.real_part() -
                        Dot(left.imaginary_part(), right.imaginary_part()),
                    left.real_part() * right.imaginary_part() +
                        right.real_part() * left.imaginary_part() +
                        Cross(left.imaginary_part(), right.imaginary_part()));
}

inline Quaternion operator/(Quaternion const& left, Quaternion const& right) {
  return left * right.Inverse();
}

inline Quaternion operator*(double const left, Quaternion const& right) {
  return Quaternion(left * right.real_part(),
                    left * right.imaginary_part());
}

inline Quaternion operator*(Quaternion const& left, double const right) {
  return Quaternion(left.real_part() * right,
                    left.imaginary_part() * right);
}

inline Quaternion operator/(Quaternion const& left, double const right) {
  return Quaternion(left.real_part() / right,
                    left.imaginary_part() / right);
}

inline std::ostream& operator<<(std::ostream& out,
                                Quaternion const& quaternion) {
  return out << quaternion.real_part() << " + "
             << quaternion.imaginary_part().x << " i + "
             << quaternion.imaginary_part().y << " j + "
             << quaternion.imaginary_part().z << " k";
}

}  // namespace geometry
}  // namespace principia
