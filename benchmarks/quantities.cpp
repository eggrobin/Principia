
// .\Release\benchmarks.exe --benchmark_filter=Cosine --benchmark_min_time=3 --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/02-13:40:44
// Benchmark                                       Time(ns)    CPU(ns) Iterations                       // NOLINT(whitespace/line_length)
// ------------------------------------------------------------------------------                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            257337     252975       2405                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            235884     246118       2472                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            249963     245522       2478                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            242904     246118       2472                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            256793     248938       2444                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform_mean       248510     247903      12271                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform_stddev       8235       2772      12271                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                      13         46   13241697                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                      21         45   13384982                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                      14         45   13387869                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                      16         45   13569181                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                      13         45   13372639                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform_mean                 16         45   66956368                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform_stddev                3          0   66956368                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                     14         46   13198184                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                     13         45   13396072                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                     15         46   13181218                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                     13         45   13385235                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                     11         45   13383427                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform_mean                13         46   66544136                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform_stddev               1          0   66544136                       // NOLINT(whitespace/line_length)

#include "benchmarks/quantities.hpp"

#include<vector>

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {
namespace benchmarks {

static void BM_DimensionfulDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)>
  std::vector<quantities::Momentum> output;
  while (state.KeepRunning()) {
    DimensionfulDiscreteCosineTransform(&output);
  }
}
BENCHMARK(BM_DimensionfulDiscreteCosineTransform);

static void BM_DoubleDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)>
  std::vector<double> output;
  while (state.KeepRunning()) {
    DoubleDiscreteCosineTransform(&output);
  }
}
BENCHMARK(BM_DoubleDiscreteCosineTransform);

static void BM_TestTypeiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)>
  std::vector<TestType> output;
  while (state.KeepRunning()) {
    TestTypeDiscreteCosineTransform(&output);
  }
}
BENCHMARK(BM_TestTypeiscreteCosineTransform);


}  // namespace benchmarks
}  // namespace principia
