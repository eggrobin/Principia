#pragma once

#include <functional>

#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/rotating_body.hpp"

namespace principia {

namespace physics {

template<typename Frame>
FrameField<Frame> CoordinateFrame() {
  return [](Position<Frame> const&) -> Rotation<Frame, Frame> {
           return Rotation<Frame, Frame>::Identity();
         };
}

}  // namespace physics
}  // namespace principia
