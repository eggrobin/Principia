﻿
#pragma once

#include <experimental/filesystem>
#include <fstream>

#include "base/not_null.hpp"
#include "serialization/journal.pb.h"

namespace principia {
namespace base {
template <typename OtherPointer> class not_null;
}  // namespace base
}  // namespace principia

namespace principia {
namespace serialization {
class Method;
}  // namespace serialization
}  // namespace principia

namespace principia {
namespace journal {

class Recorder {
 public:
  explicit Recorder(std::experimental::filesystem::path const& path);
  ~Recorder();

  void Write(serialization::Method const& method);

  static void Activate(base::not_null<Recorder*> const recorder);
  static void Deactivate();
  static bool IsActivated();

 private:
  std::ofstream stream_;

  static Recorder* active_recorder_;

  template<typename>
  friend class Method;
  friend class RecorderTest;
};

}  // namespace journal
}  // namespace principia
