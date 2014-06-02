
// .\Release\benchmarks.exe --benchmark_filter=Cosine --benchmark_min_time=3 --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/02-16:48:10
// Benchmark                                       Time(ns)    CPU(ns) Iterations                       // NOLINT(whitespace/line_length)
// ------------------------------------------------------------------------------                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            256409     252555       2409                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            235025     245621       2477                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            249719     245324       2480                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            242947     246118       2472                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            256617     249244       2441                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform_mean       248075     247742      12279                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform_stddev       8266       2759      12279                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     488        520    1169638                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     458        512    1188848                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     487        511    1190203                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     502        518    1174661                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     458        512    1187827                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform_mean                478        515    5911177                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform_stddev               17          4    5911177                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    496        526    1155744                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    491        517    1176484                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    517        541    1125406                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    497        521    1168319                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    473        507    1200313                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform_mean               495        522    5826266                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform_stddev              14         11    5826266                       // NOLINT(whitespace/line_length)

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
