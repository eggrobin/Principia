#pragma once

#include "base/macros.hpp"
#include "base/not_null.hpp"

#if HAS_SHARED_MUTEX
#include <shared_mutex>
#else
#include <mutex>
#endif

namespace principia {
namespace base {

#if HAS_SHARED_MUTEX

using shared_mutex = std::shared_mutex;

// A helper class that the language designer didn't think useful of providing.
template<typename Mutex>
class SCOPED_CAPABILITY shared_lock_guard final {
 public:
  explicit shared_lock_guard(Mutex& mutex) ACQUIRE_SHARED(mutex);
  ~shared_lock_guard() RELEASE();

 private:
  not_null<Mutex*> const mutex_;
};

#else

using shared_mutex = std::mutex;
template<typename Mutex>
using shared_lock_guard = std::lock_guard<Mutex>;

#endif

}  // namespace base
}  // namespace principia

#include "base/shared_lock_guard_body.hpp"
