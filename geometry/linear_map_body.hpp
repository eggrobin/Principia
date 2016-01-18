﻿#pragma once

#include "base/not_null.hpp"
#include "geometry/linear_map.hpp"
#include "serialization/geometry.pb.h"

namespace principia {
namespace geometry {

template<typename FromFrame, typename ToFrame>
void LinearMap<FromFrame, ToFrame>::WriteToMessage(
    not_null<serialization::LinearMap*> const message) {
  FromFrame::WriteToMessage(message->mutable_from_frame());
  ToFrame::WriteToMessage(message->mutable_to_frame());
}

template<typename FromFrame, typename ToFrame>
void LinearMap<FromFrame, ToFrame>::ReadFromMessage(
    serialization::LinearMap const& message) {
  FromFrame::ReadFromMessage(message.from_frame());
  ToFrame::ReadFromMessage(message.to_frame());
}

}  // namespace geometry
}  // namespace principia
