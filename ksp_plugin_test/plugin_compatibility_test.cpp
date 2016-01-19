﻿
#include <stdint.h>
#include <experimental/filesystem>
#include <fstream>  // IWYU pragma: keep
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "google/protobuf/io/coded_stream.h"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/plugin.hpp"
#include "physics/barycentric_rotating_dynamic_frame.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/oblate_body.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using base::Array;
using base::UniqueBytes;
using base::HexadecimalDecode;
using geometry::Bivector;
using geometry::Trivector;
using geometry::Vector;
using quantities::Length;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

namespace ksp_plugin {

class TestablePlugin : public Plugin {
 public:
  void KeepAllVessels();

  static not_null<std::unique_ptr<TestablePlugin>> ReadFromMessage(
    serialization::Plugin const& message);
};

void TestablePlugin::KeepAllVessels() {
  for (auto const& pair : vessels_) {
    auto const& vessel = pair.second;
    kept_vessels_.insert(vessel.get());
  }
}

not_null<std::unique_ptr<TestablePlugin>> TestablePlugin::ReadFromMessage(
  serialization::Plugin const& message) {
  std::unique_ptr<Plugin> plugin = Plugin::ReadFromMessage(message);
  return std::unique_ptr<TestablePlugin>(
      reinterpret_cast<TestablePlugin*>(plugin.release()));
}

class PluginCompatibilityTest : public testing::Test {
};

TEST_F(PluginCompatibilityTest, PreBorel) {
  serialization::Multivector message;

  Vector<Length, Barycentric> const v({ -1 * Metre, 2 * Metre, 3 * Metre });
  v.WriteToMessage(&message);
  message.mutable_frame()->set_tag(serialization::Frame::PRE_BOREL_BARYCENTRIC);
  Vector<Length, Barycentric> const w =
    Vector<Length, Barycentric>::ReadFromMessage(message);
  Vector<Length, Barycentric> const expected_w(
  { -1 * Metre, 3 * Metre, 2 * Metre });
  EXPECT_EQ(expected_w, w);

  Bivector<Length, Barycentric> const b({ 4 * Metre, 5 * Metre, -6 * Metre });
  b.WriteToMessage(&message);
  message.mutable_frame()->set_tag(serialization::Frame::PRE_BOREL_BARYCENTRIC);
  Bivector<Length, Barycentric> const c =
    Bivector<Length, Barycentric>::ReadFromMessage(message);
  Bivector<Length, Barycentric> const expected_c(
  { -4 * Metre, 6 * Metre, -5 * Metre });
  EXPECT_EQ(expected_c, c);

  Trivector<Length, Barycentric> const t(-7 * Metre);
  t.WriteToMessage(&message);
  message.mutable_frame()->set_tag(serialization::Frame::PRE_BOREL_BARYCENTRIC);
  Trivector<Length, Barycentric> const u =
    Trivector<Length, Barycentric>::ReadFromMessage(message);
  Trivector<Length, Barycentric> const expected_u(7 * Metre);
  EXPECT_EQ(expected_u, u);
}

TEST_F(PluginCompatibilityTest, PreBourbaki) {
  // Read the entire hex data.
  std::fstream file = std::fstream(
      SOLUTION_DIR / "ksp_plugin_test" / "pre_bourbaki.proto.hex");
  CHECK(file.good());
  std::string hex;
  while (!file.eof()) {
    std::string line;
    std::getline(file, line);
    for (auto const c : line) {
      if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')) {
        hex.append(1, c);
      }
    }
  }
  file.close();

  // Parse it and convert to binary data.
  UniqueBytes bin(hex.size() / 2);
  HexadecimalDecode(Array<std::uint8_t const>(
                        reinterpret_cast<const std::uint8_t*>(hex.c_str()),
                        hex.size()),
                    bin.get());

  // Construct a protocol buffer from the binary data.
  google::protobuf::io::CodedInputStream coded_input_stream(
      bin.data.get(), static_cast<int>(bin.size));
  serialization::Plugin pre_bourbaki_serialized_plugin;
  CHECK(pre_bourbaki_serialized_plugin.MergeFromCodedStream(
            &coded_input_stream));

  // Construct a plugin from the protocol buffer.
  auto plugin = TestablePlugin::ReadFromMessage(pre_bourbaki_serialized_plugin);

  // Do some operations on the plugin.
  plugin->KeepAllVessels();
  plugin->AdvanceTime(plugin->CurrentTime() + 1 * Second, 2 * Radian);
  plugin->AdvanceTime(plugin->CurrentTime() + 1 * Hour, 3 * Radian);

  // Serialize and deserialize it in the new format.
  serialization::Plugin post_bourbaki_serialized_plugin;
  plugin->WriteToMessage(&post_bourbaki_serialized_plugin);
  plugin = TestablePlugin::ReadFromMessage(post_bourbaki_serialized_plugin);
}

}  // namespace ksp_plugin
}  // namespace principia
