
#pragma once

#include "base/sink_source.hpp"

#include <algorithm>

namespace principia {
namespace base {
namespace internal_sink_source {

template<typename Element>
principia::base::internal_sink_source::ArraySource<Element>::ArraySource(
    Array<Element> const& array)
    : array_(array) {}

template<typename Element>
size_t principia::base::internal_sink_source::ArraySource<Element>::Available()
    const {
  return array_.size - next_to_read_;
}

template<typename Element>
const char* ArraySource<Element>::Peek(size_t* const length) {
  *length = array_.size - next_to_read_;
  return reinterpret_cast<const char*>(array_.data + next_to_read_);
}

template<typename Element>
void ArraySource<Element>::Skip(size_t const n) {
  next_to_read_ += n;
}

template<typename Element>
principia::base::internal_sink_source::ArraySink<Element>::ArraySink(
    Array<Element> const& array)
    : array_(array) {}

template<typename Element>
Array<Element> ArraySink<Element>::array() const {
  Array<Element> result;
  result.data = array_.data;
  result.size = next_to_write_;
  return result;
}

template<typename Element>
void principia::base::internal_sink_source::ArraySink<Element>::Append(
    const char* const data,
    size_t const n) {
  // Do no copying if the caller filled in the result of GetAppendBuffer()
  if (data != reinterpret_cast<const char*>(array_.data + next_to_write_)) {
    memcpy(array_.data + next_to_write_, data, n);
  }
  next_to_write_ += n;
}

template<typename Element>
char* ArraySink<Element>::GetAppendBuffer(
    size_t const min_size,
    size_t const desired_size_hint,
    char* const scratch,
    size_t const scratch_size,
    size_t* const allocated_size) {
  *allocated_size = std::min(static_cast<std::int64_t>(desired_size_hint),
                             array_.size - next_to_write_);
  return reinterpret_cast<char*>(array_.data + next_to_write_);
}

}  // namespace internal_sink_source
}  // namespace base
}  // namespace principia
