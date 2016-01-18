﻿
#pragma once

#include <ostream>
#include <string>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "geometry/sign.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

class Sign;

template<typename Scalar>
Sign::Sign(Scalar const& scalar) : negative_(scalar < Scalar{}) {}

inline bool Sign::Negative() const {
  return negative_;
}

inline bool Sign::Positive() const {
  return !negative_;
}

inline bool Sign::operator==(Sign const& other) const {
  return negative_ == other.negative_;
}

inline bool Sign::operator!=(Sign const& other) const {
  return negative_ != other.negative_;
}

inline void Sign::WriteToMessage(
    not_null<serialization::Sign*> const message) const {
  message->set_negative(negative_);
}

inline Sign Sign::ReadFromMessage(serialization::Sign const& message) {
  return Sign(message.negative() ? -1 : 1);
}

inline Sign operator*(Sign const& left, Sign const& right) {
  return Sign(left.negative_ == right.negative_ ? 1 : -1);
}

template<typename T>
FORCE_INLINE T operator*(Sign const& left, T const& right) {
  return left.negative_ ? -right : right;
}

inline std::string DebugString(Sign const& sign) {
  return sign.Negative() ? "-" : "+";
}

inline std::ostream& operator<<(std::ostream& out, Sign const& sign) {
  out << DebugString(sign);
  return out;
}

}  // namespace geometry
}  // namespace principia
