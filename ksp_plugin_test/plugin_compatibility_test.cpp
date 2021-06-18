
#include <filesystem>
#include <fstream>
#include <map>
#include <string>

#include "base/array.hpp"
#include "base/hexadecimal.hpp"
#include "geometry/grassmann.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "google/protobuf/io/coded_stream.h"
#include "gtest/gtest.h"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/plugin.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/serialization.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_plugin {

using base::Array;
using base::HexadecimalEncoder;
using base::UniqueArray;
using geometry::Bivector;
using geometry::Trivector;
using geometry::Vector;
using quantities::Length;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using ::testing::AllOf;
using ::testing::AnyOf;
using ::testing::Gt;
using ::testing::Lt;

class TestablePlugin : public Plugin {
 public:
  void KeepAllVessels();

  std::map<GUID, not_null<Vessel const*>> vessels() const;

  static not_null<std::unique_ptr<TestablePlugin>> ReadFromMessage(
      serialization::Plugin const& message);
};

void TestablePlugin::KeepAllVessels() {
  for (auto const& [_, vessel] : vessels_) {
    kept_vessels_.insert(vessel.get());
  }
}

std::map<GUID, not_null<Vessel const*>> TestablePlugin::vessels() const {
  std::map<GUID, not_null<Vessel const*>> result;
  for (auto const& [guid, vessel] : vessels_) {
    result.insert(std::make_pair(guid, vessel.get()));
  }
  return result;
}

not_null<std::unique_ptr<TestablePlugin>> TestablePlugin::ReadFromMessage(
    serialization::Plugin const& message) {
  std::unique_ptr<Plugin> plugin = Plugin::ReadFromMessage(message);
  return std::unique_ptr<TestablePlugin>(
      static_cast<TestablePlugin*>(plugin.release()));
}

class PluginCompatibilityTest : public testing::Test {
 protected:
  serialization::Plugin ReadFromFile(std::string const& filename) {
    // Open the file and read hexadecimal data.
    std::fstream file =
        std::fstream(SOLUTION_DIR / "ksp_plugin_test" / filename);
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

    // Parse the hexadecimal data and convert it to binary data.
    HexadecimalEncoder</*null_terminated=*/false> encoder;
    auto const bin = encoder.Decode(hex);

    // Construct a protocol buffer from the binary data.
    google::protobuf::io::CodedInputStream coded_input_stream(
        bin.data.get(), static_cast<int>(bin.size));
    serialization::Plugin message;
    CHECK(message.MergeFromCodedStream(&coded_input_stream));

    return message;
  }
};

TEST_F(PluginCompatibilityTest, PreCartan) {
  // This space for rent.
}

TEST_F(PluginCompatibilityTest, Reach) {
  std::unique_ptr<Plugin const> const plugin = []() {
    LOG(ERROR) << u8"Reading file…";
    auto const lines = testing_utilities::ReadLinesFromBase64File(
        SOLUTION_DIR / "ksp_plugin_test" / "Reach.sfs");
    base::PushDeserializer* deserializer = nullptr;
    Plugin const* plugin;
    LOG(ERROR) << lines.size() << u8" lines. Deserializing…";
    for (std::string const& line : lines) {
      interface::principia__DeserializePlugin(line.c_str(),
                                              &deserializer,
                                              &plugin,
                                              /*compressor=*/"gipfeli",
                                              "base64");
    }
    interface::principia__DeserializePlugin("",
                                            &deserializer,
                                            &plugin,
                                            /*compressor=*/"gipfeli",
                                            "base64");
    LOG(ERROR) << "Deserialization complete.";
    return std::unique_ptr<Plugin const>(plugin);
  }();
  auto const test =
      plugin->GetVessel("f2d77873-4776-4809-9dfb-de9e7a0620a6");  // TEST.
  auto const infnity =
      plugin->GetVessel("29142a79-7acd-47a9-a34d-f9f2a8e1b4ed");  // IFNITY-5.2.
  for (Vessel const* vessel : {test, infnity}) {
    LOG(ERROR) << vessel->name() << ":";
    if (vessel->has_flight_plan()) {
      auto const& flight_plan = vessel->flight_plan();
      LOG(ERROR) << flight_plan.number_of_manœuvres() << " manœuvres";
      LOG(ERROR) << "Flight plan initial time: " << flight_plan.initial_time();
      for (int i = 0; i < flight_plan.number_of_manœuvres(); ++i) {
        LOG(ERROR) << flight_plan.GetManœuvre(i).Δv().Norm() << " at "
                   << flight_plan.GetManœuvre(i).initial_time();
      }
    } else {
      LOG(ERROR) << "No flight plan.";
    }
    LOG(ERROR) << "Psychohistory range: "
               << vessel->psychohistory().front().time << " .. "
               << vessel->psychohistory().back().time;
  }
}

}  // namespace internal_plugin
}  // namespace ksp_plugin
}  // namespace principia
