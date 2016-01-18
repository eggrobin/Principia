﻿#pragma once

#include <memory>
#include <ostream>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "glog/logging.h"
#include "physics/body.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {

class Body;
template <typename Frame> class OblateBody;

template<typename Frame>
bool Body::is_compatible_with() const {
  return CompatibilityHelper<Frame,
                             Frame::is_inertial>::is_compatible_with(this);
}

inline not_null<std::unique_ptr<Body>> Body::ReadFromMessage(
    serialization::Body const& message) {
  if (message.has_massless_body()) {
    return MasslessBody::ReadFromMessage(message.massless_body());
  } else if (message.has_massive_body()) {
    return MassiveBody::ReadFromMessage(message.massive_body());
  } else {
    LOG(FATAL) << "Body is neither massive nor massless";
    base::noreturn();
  }
}

template<typename Frame>
class Body::CompatibilityHelper<Frame, false> {
 public:
  static bool is_compatible_with(not_null<Body const*> const body);
};

template<typename Frame>
class Body::CompatibilityHelper<Frame, true> {
 public:
  static bool is_compatible_with(not_null<Body const*> const body);
};

template<typename Frame>
bool Body::CompatibilityHelper<Frame, false>::is_compatible_with(
    not_null<Body const*> const body) {
  return !body->is_oblate();
}

template<typename Frame>
bool Body::CompatibilityHelper<Frame, true>::is_compatible_with(
    not_null<Body const*> const body) {
  return !body->is_oblate() ||
         dynamic_cast<OblateBody<Frame> const*>(
            static_cast<Body const*>(body)) != nullptr;
}

}  // namespace physics
}  // namespace principia
