﻿
#include <string.h>
#include <string>

#include "base/array.hpp"
#include "base/get_line.hpp"
#include "base/hexadecimal.hpp"
#include "experimental/filesystem"
#include "glog/logging.h"
#include "journal/player.hpp"
#include "serialization/journal.pb.h"

namespace principia {

using base::GetLine;
using base::HexadecimalDecode;
using base::UniqueBytes;

namespace journal {

Player::Player(std::experimental::filesystem::path const& path)
    : stream_(path, std::ios::in) {
  CHECK(!stream_.fail());
}

bool Player::Play() {
  std::unique_ptr<serialization::Method> method = Read();
  if (method == nullptr) {
    return false;
  }

#include "journal/player.generated.cc"

  return true;
}

std::unique_ptr<serialization::Method> Player::Read() {
  std::string const line = GetLine(&stream_);
  if (line.empty()) {
    return nullptr;
  }

  uint8_t const* const hexadecimal =
      reinterpret_cast<uint8_t const*>(line.c_str());
  int const hexadecimal_size = strlen(line.c_str());
  UniqueBytes bytes(hexadecimal_size >> 1);
  HexadecimalDecode({hexadecimal, hexadecimal_size},
                    {bytes.data.get(), bytes.size});
  auto method = std::make_unique<serialization::Method>();
  CHECK(method->ParseFromArray(bytes.data.get(),
                               static_cast<int>(bytes.size)));

  return method;
}

}  // namespace journal
}  // namespace principia
