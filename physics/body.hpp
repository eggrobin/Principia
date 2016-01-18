﻿
#ifndef PRINCIPIA_PHYSICS_BODY_HPP_
#define PRINCIPIA_PHYSICS_BODY_HPP_

#include <memory>

#include "base/not_null.hpp"
#include "serialization/physics.pb.h"

namespace principia {

namespace serialization {
class Body;
}  // namespace serialization

using base::not_null;

namespace physics {

class Body {
 public:
  ~Body() = default;

  // Returns true iff this body is massless.
  virtual bool is_massless() const = 0;

  // Returns true iff this body is oblate (which implies massive).
  virtual bool is_oblate() const = 0;

  // Returns true iff this body is compatible with the given frame (either
  // because it is spherical or because its axis is expressed in the same
  // frame).
  template<typename Frame>
  bool is_compatible_with() const;

  virtual void WriteToMessage(not_null<serialization::Body*> message) const = 0;

  // Dispatches to one of the subclasses depending on the contents of the
  // message.
  // Beware the Jabberwock, my son!  If it dispatches to |RotatingBody|, this
  // method will return an |RotatingBody<UnknownFrame>|.  Use |reinterpret_cast|
  // afterwards as appropriate if the frame is known.
  static not_null<std::unique_ptr<Body>> ReadFromMessage(
      serialization::Body const& message);

 protected:
  Body() = default;

 private:
  // A helper class which is here just so that we can specialize it on
  // |is_inertial|.  This is necessary because we cannot dynamic cast to
  // OblateBody<Frame> if |Frame| is not inertial.
  template<typename Frame, bool is_inertial>
  class CompatibilityHelper {
   public:
     static bool is_compatible_with(not_null<Body const*> const body);
   private:
     CompatibilityHelper() = delete;
  };
};

}  // namespace physics
}  // namespace principia

#include "physics/body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BODY_HPP_
