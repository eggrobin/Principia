﻿
#pragma once

// This code comes from:
// https://github.com/google/protobuf/tree/master/src/google/protobuf/stubs
// and was adapted to Visual Studio and to the needs of this project.

// Protocol Buffers - Google's data interchange format
// Copyright 2008 Google Inc.  All rights reserved.
// https://developers.google.com/protocol-buffers/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iosfwd>
#include <string>

#include "glog/logging.h"

// TODO(phl): Many of the functions in this file should be made constexpr.
// Also, we should use string_view.

namespace principia {
namespace base {

// See https://cloud.google.com/vision/reference/rest/v1/Code for recommended
// usage of these codes.
enum class Error {
  OK = 0,
  CANCELLED = 1,
  UNKNOWN = 2,
  INVALID_ARGUMENT = 3,
  DEADLINE_EXCEEDED = 4,
  NOT_FOUND = 5,
  ALREADY_EXISTS = 6,
  PERMISSION_DENIED = 7,
  UNAUTHENTICATED = 16,
  RESOURCE_EXHAUSTED = 8,
  FAILED_PRECONDITION = 9,
  ABORTED = 10,
  OUT_OF_RANGE = 11,
  UNIMPLEMENTED = 12,
  INTERNAL = 13,
  UNAVAILABLE = 14,
  DATA_LOSS = 15,
};

Error operator|(Error left, Error right);
Error& operator|=(Error& left, Error right);

std::string ErrorToString(Error error);

class Status final {
 public:
  // Creates a "successful" status.
  Status() = default;

  Status(Error error, std::string const& message);

  // Some pre-defined Status objects.
  static const Status OK;
  static const Status CANCELLED;
  static const Status UNKNOWN;

  // Accessors.
  bool ok() const;
  Error error() const;
  std::string const& message() const;

  bool operator==(Status const& s) const;
  bool operator!=(Status const& s) const;

  void Update(Status const& s);

  // Returns a combination of the error code name and message.
  std::string ToString() const;

 private:
  Error error_ = Error::OK;
  std::string message_;
};

inline Status const& GetStatus(Status const& s) {
  return s;
}

// Prints a human-readable representation of |s| to |os|.
std::ostream& operator<<(std::ostream& os, Status const& s);

#define CHECK_OK(value) CHECK_EQ((value), ::principia::base::Status::OK)
#define EXPECT_OK(value) \
  EXPECT_THAT((value), ::principia::testing_utilities::IsOk());

#define RETURN_IF_ERROR(expr)                                                \
  do {                                                                       \
    /* Using _status below to avoid capture problems if expr is "status". */ \
    ::principia::base::Status const _status =                                \
        (::principia::base::GetStatus(expr));                                \
    if (!_status.ok())                                                       \
      return _status;                                                        \
  } while (false)

}  // namespace base
}  // namespace principia
