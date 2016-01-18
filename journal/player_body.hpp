﻿#pragma once

#include <list>

#include "google/protobuf/extension_set.h"
#include "journal/player.hpp"

namespace principia {
namespace serialization {
class Method;
}  // namespace serialization
}  // namespace principia

namespace principia {
namespace journal {

template<typename Profile>
bool Player::RunIfAppropriate(serialization::Method const& method) {
  if (method.HasExtension(Profile::Message::extension)) {
    Profile::Run(method.GetExtension(Profile::Message::extension),
                 &pointer_map_);
    return true;
  }
  return false;
}

}  // namespace journal
}  // namespace principia
