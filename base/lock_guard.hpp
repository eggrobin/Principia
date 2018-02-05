#pragma once

#include <mutex>

#include "base/macros.hpp"

namespace principia {
namespace base {

// A drop-in replacement for (and a trivial wrapper around) |std::lock_guard|,
// with thread annotations.
template<typename Mutex>
class SCOPED_CAPABILITY lock_guard {
 public:
  explicit lock_guard(Mutex& mutex) ACQUIRE(mutex) : guard_(mutex) {}
  ~lock_guard() RELEASE() = default;
 private:
  std::lock_guard<Mutex> guard_;
};

// A |unique_lock| which really is just a |lock_guard|, except that it can be
// used to |wait| on a |std::condition_variable|.
// It cannot be unlocked and relocked (thread safety analysis does not support
// reacquisition of scoped capabilities), and it cannot be moved.
// See http://lists.llvm.org/pipermail/cfe-dev/2016-November/051468.html.
template<typename Mutex>
class SCOPED_CAPABILITY unique_lock {
 public:
  explicit unique_lock(Mutex& mutex) ACQUIRE(mutex) : lock_(mutex) {}
  ~unique_lock() RELEASE() = default;

  // No support for moving locks in thread safety analysis.
  unique_lock(unique_lock&&) = delete;
  unique_lock& operator=(unique_lock&&) = delete;

  // Calling the predicate ourselves allows the |REQUIRES| on the Predicate to
  // be checked by thread safety analysis.
  template<typename Predicate>
  void wait(std::condition_variable& condition_variable, Predicate predicate) {
    while (!predicate()) {
      condition_variable.wait(lock_);
    }
  }

 private:
  std::unique_lock<Mutex> lock_;
};

}  // namespace base
}  // namespace principia
