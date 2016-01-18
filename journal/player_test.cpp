﻿
#include <functional>
#include <sstream>
#include <string>

#include "benchmark/benchmark_api.h"
#include "experimental/filesystem"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "journal/method.hpp"
#include "journal/player.hpp"
#include "journal/profiles.hpp"
#include "journal/recorder.hpp"
#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {
namespace journal {

void BM_PlayForReal(benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    Player player("P:\\Public Mockingbird\\Principia\\JOURNAL.20151206-170008");
    int count = 0;
    while (player.Play()) {
      ++count;
      LOG_IF(ERROR, (count % 100'000) == 0)
          << count << " journal entries replayed";
    }
  }
}

#if 0
BENCHMARK(BM_PlayForReal);
#endif

class PlayerTest : public testing::Test {
 protected:
  static void SetUpTestCase() {
    benchmark::RunSpecifiedBenchmarks();
  }

  PlayerTest()
      : test_name_(
            testing::UnitTest::GetInstance()->current_test_info()->name()),
        plugin_(interface::principia__NewPlugin(1, 2)),
        recorder_(new Recorder(test_name_ + ".journal.hex")) {
    Recorder::Activate(recorder_);
  }

  ~PlayerTest() override {
    Recorder::Deactivate();
  }

  std::string const test_name_;
  std::unique_ptr<ksp_plugin::Plugin> plugin_;
  Recorder* recorder_;
};

TEST_F(PlayerTest, PlayTiny) {
  {
    Method<NewPlugin> m({1, 2});
    m.Return(plugin_.get());
  }
  {
    const ksp_plugin::Plugin* plugin = plugin_.get();
    Method<DeletePlugin> m({&plugin}, {&plugin});
    m.Return();
  }

  Player player(test_name_ + ".journal.hex");

  // Replay the journal.  Note that the journal doesn't grow as we replay
  // because we didn't call principia__ActivateRecorder so there is no active
  // journal in the ksp_plugin assembly.
  int count = 0;
  while (player.Play()) {
    ++count;
  }
  EXPECT_EQ(2, count);
}

}  // namespace journal
}  // namespace principia
