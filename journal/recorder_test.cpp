﻿#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "experimental/filesystem"
#include "google/protobuf/extension_set.h"
#include "gtest/gtest-death-test.h"
#include "gtest/gtest.h"
#include "journal/method.hpp"
#include "journal/player.hpp"
#include "journal/profiles.hpp"
#include "journal/recorder.hpp"
#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/plugin.hpp"
#include "serialization/journal.pb.h"

namespace principia {
namespace journal {

class RecorderTest : public testing::Test {
 protected:
  RecorderTest()
      : test_name_(
            testing::UnitTest::GetInstance()->current_test_info()->name()),
        plugin_(interface::principia__NewPlugin(1, 2)),
        recorder_(new Recorder(test_name_ + ".journal.hex")) {
    Recorder::Activate(recorder_);
  }

  ~RecorderTest() override {
    Recorder::Deactivate();
  }

  static std::vector<serialization::Method> ReadAll(
      std::experimental::filesystem::path const& path) {
    std::vector<serialization::Method> methods;
    Player player(path);
    for (std::unique_ptr<serialization::Method> method = player.Read();
         method != nullptr;
         method = player.Read()) {
      methods.push_back(*method);
    }
    return methods;
  }


  std::string const test_name_;
  std::unique_ptr<ksp_plugin::Plugin> plugin_;
  Recorder* recorder_;
};

using JournalDeathTest = RecorderTest;

TEST_F(JournalDeathTest, Return) {
  EXPECT_DEATH({
    Method<NewPlugin> m({1, 2});
    m.Return(plugin_.get());
    m.Return();
  },
  "!returned_");
  EXPECT_DEATH({
    const ksp_plugin::Plugin* plugin = plugin_.get();
    Method<DeletePlugin> m({&plugin}, {&plugin});
    m.Return();
    m.Return();
  },
  "!returned_");
  EXPECT_DEATH({
    Method<NewPlugin> m({1, 2});
  },
  "returned_");
}

TEST_F(RecorderTest, Recording) {
  {
    const ksp_plugin::Plugin* plugin = plugin_.get();
    Method<DeletePlugin> m({&plugin}, {&plugin});
    m.Return();
  }
  {
    Method<NewPlugin> m({1, 2});
    m.Return(plugin_.get());
  }

  std::vector<serialization::Method> const methods =
      ReadAll(test_name_ + ".journal.hex");
  EXPECT_EQ(2, methods.size());
  auto it = methods.begin();
  {
    EXPECT_TRUE(it->HasExtension(serialization::DeletePlugin::extension));
    auto const& extension =
        it->GetExtension(serialization::DeletePlugin::extension);
    EXPECT_TRUE(extension.has_in());
    EXPECT_NE(0, extension.in().plugin());
    EXPECT_TRUE(extension.has_out());
    EXPECT_NE(0, extension.out().plugin());
    EXPECT_EQ(extension.in().plugin(), extension.out().plugin());
  }
  ++it;
  {
    EXPECT_TRUE(it->HasExtension(serialization::NewPlugin::extension));
    auto const& extension =
        it->GetExtension(serialization::NewPlugin::extension);
    EXPECT_TRUE(extension.has_in());
    EXPECT_EQ(1, extension.in().initial_time());
    EXPECT_EQ(2, extension.in().planetarium_rotation_in_degrees());
    EXPECT_TRUE(extension.has_return_());
    EXPECT_NE(0, extension.return_().result());
  }
}

}  // namespace journal
}  // namespace principia
