
#pragma once

#include "base/not_null.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace testing_utilities {

template<typename T>
not_null<T> make_not_null() {
  return reinterpret_cast<T>(0xDEADBEEF);
}

}  // namespace testing_utilities
}  // namespace principia
