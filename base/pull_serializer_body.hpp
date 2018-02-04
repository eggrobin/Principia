﻿
#pragma once

#include "base/lock_guard.hpp"
#include "base/pull_serializer.hpp"

#include <algorithm>

namespace principia {
namespace base {
namespace internal_pull_serializer {

using std::placeholders::_1;
using std::swap;

inline DelegatingArrayOutputStream::DelegatingArrayOutputStream(
    Bytes const bytes,
    std::function<Bytes(Bytes const bytes)> on_full)
    : bytes_(bytes),
      on_full_(std::move(on_full)),
      byte_count_(0),
      position_(0),
      last_returned_size_(0) {}

inline bool DelegatingArrayOutputStream::Next(void** const data,
                                              int* const size) {
  if (position_ == bytes_.size) {
    // We're at the end of the array.  Hand the current array over to the
    // callback and start filling the next one.
    bytes_ = on_full_(bytes_);
    position_ = 0;
  }
  CHECK_LT(position_, bytes_.size);
  last_returned_size_ = bytes_.size - position_;
  *data = &bytes_.data[position_];
  *size = static_cast<int>(last_returned_size_);
  byte_count_ += last_returned_size_;
  position_ += last_returned_size_;
  return true;
}

inline void DelegatingArrayOutputStream::BackUp(int count) {
  CHECK_GT(last_returned_size_, 0)
      << "BackUp() can only be called after a successful Next().";
  CHECK_LE(count, last_returned_size_);
  CHECK_GE(count, 0);
  byte_count_ -= count;
  position_ -= count;
  // This is called at the end of the stream, in which case we must notify the
  // client about any data remaining in the stream.  If this is called at other
  // times, well, notifying the client doesn't hurt as long as we don't pass a
  // size of 0.
  if (position_ > 0) {
    bytes_ = on_full_(Bytes(bytes_.data, position_));
    position_ = 0;
  }
  last_returned_size_ = 0;
}

inline std::int64_t DelegatingArrayOutputStream::ByteCount() const {
  return byte_count_;
}

inline PullSerializer::PullSerializer(int const chunk_size,
                                      int const number_of_chunks)
    : chunk_size_(chunk_size),
      number_of_chunks_(number_of_chunks),
      data_(std::make_unique<std::uint8_t[]>(chunk_size_ * number_of_chunks_)),
      stream_(Bytes(data_.get(), chunk_size_),
              std::bind(&PullSerializer::Push, this, _1)) {
  // Mark all the chunks as free except the last one which is a sentinel for the
  // |queue_|.  The 0th chunk has been passed to the stream, but it's still free
  // until the first call to |on_full|.
  for (int i = 0; i < number_of_chunks_ - 1; ++i) {
    free_.push(data_.get() + i * chunk_size_);
  }
  queue_.push(Bytes(data_.get() + (number_of_chunks_ - 1) * chunk_size_, 0));
}

inline PullSerializer::~PullSerializer() {
  if (thread_ != nullptr) {
    thread_->join();
  }
}

inline void PullSerializer::Start(
    not_null<std::unique_ptr<google::protobuf::Message const>> message) {
  CHECK(thread_ == nullptr);
  message_ = std::move(message);
  thread_ = std::make_unique<std::thread>([this](){
    CHECK(message_->SerializeToZeroCopyStream(&stream_));
    // Put a sentinel at the end of the serialized stream so that the client
    // knows that this is the end.
    Bytes bytes;
    {
      base::unique_lock<std::mutex> l(lock_);
      CHECK(!free_.empty());
      bytes = Bytes(free_.front(), 0);
    }
    Push(bytes);
  });
}

inline Bytes PullSerializer::Pull() {
  Bytes result;
  {
    base::unique_lock<std::mutex> l(lock_);
    // The element at the front of the queue is the one that was last returned
    // by |Pull| and must be dropped and freed.
    queue_has_elements_.wait(l, [this]() { return queue_.size() > 1; });
    CHECK_LE(2u, queue_.size());
    free_.push(queue_.front().data);
    queue_.pop();
    result = queue_.front();
    CHECK_EQ(number_of_chunks_, queue_.size() + free_.size());
  }
  queue_has_room_.notify_all();
  return result;
}

inline Bytes PullSerializer::Push(Bytes const bytes) {
  Bytes result;
  CHECK_GE(chunk_size_, bytes.size);
  {
    base::unique_lock<std::mutex> l(lock_);
    queue_has_room_.wait(l, [this]() {
      // -1 here is because we want to ensure that there is an entry in the
      // (real) free list.
      return queue_.size() < static_cast<std::size_t>(number_of_chunks_) - 1;
    });
    queue_.emplace(bytes.data, bytes.size);
    CHECK_LE(2u, free_.size());
    CHECK_EQ(free_.front(), bytes.data);
    free_.pop();
    result = Bytes(free_.front(), chunk_size_);
    CHECK_EQ(number_of_chunks_, queue_.size() + free_.size());
  }
  queue_has_elements_.notify_all();
  return result;
}

}  // namespace internal_pull_serializer
}  // namespace base
}  // namespace principia
