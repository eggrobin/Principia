
#pragma once

#include "base/lock_guard.hpp"
#include "base/thread_pool.hpp"

namespace principia {
namespace base {
namespace internal_thread_pool {

// A helper function that is specialized for void because void is not really a
// type.
template<typename T>
void ExecuteAndSetValue(std::function<T()> const& function,
                        std::promise<T>& promise) {
  promise.set_value(function());
}

template<>
inline void ExecuteAndSetValue<void>(std::function<void()> const& function,
                                     std::promise<void>& promise) {
  function();
  promise.set_value();
}

}  // namespace internal_thread_pool

template<typename T>
ThreadPool<T>::ThreadPool(std::int64_t const pool_size) {
  for (std::int64_t i = 0; i < pool_size; ++i) {
    threads_.emplace_back(std::bind(&ThreadPool::DequeueCallAndExecute, this));
  }
}

template<typename T>
ThreadPool<T>::~ThreadPool() {
  {
    base::lock_guard<std::mutex> l(lock_);
    shutdown_ = true;
  }
  has_calls_or_shutdown_.notify_all();
  for (auto& thread : threads_) {
    thread.join();
  }
}

template<typename T>
std::future<T> ThreadPool<T>::Add(std::function<T()> function) {
  std::future<T> result;
  {
    base::lock_guard<std::mutex> l(lock_);
    calls_.push_back({std::move(function), std::promise<T>()});
    result = calls_.back().promise.get_future();
  }
  has_calls_or_shutdown_.notify_one();
  return result;
}

template<typename T>
void ThreadPool<T>::DequeueCallAndExecute() {
  for (;;) {
    Call this_call;

    // Wait until either the queue contains an element or this class is shutting
    // down.
    {
      base::unique_lock<std::mutex> l(lock_);
      has_calls_or_shutdown_.wait(l,
                                  [this] {
                                      return shutdown_ || !calls_.empty();
                                  });
      if (shutdown_) {
        break;
      }
      this_call = std::move(calls_.front());
      calls_.pop_front();
    }

    // Execute the function without holding the |lock_| as it might take some
    // time.
    internal_thread_pool::ExecuteAndSetValue(this_call.function,
                                             this_call.promise);
  }
}

}  // namespace base
}  // namespace principia
