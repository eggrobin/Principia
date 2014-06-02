
// .\Release\benchmarks.exe --benchmark_filter=Cosine --benchmark_min_time=3 --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/02-16:56:08
// Benchmark                                       Time(ns)    CPU(ns) Iterations                       // NOLINT(whitespace/line_length)
// ------------------------------------------------------------------------------                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            285523     293065       2076                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            301929     293207       2075                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            294263     289440       2102                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            277445     288891       2106                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform            294944     290131       2097                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform_mean       290798     290935      10456                       // NOLINT(whitespace/line_length)
// BM_DimensionfulDiscreteCosineTransform_stddev       8479       1829      10456                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     537        569    1068981                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     526        549    1108009                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     525        548    1110685                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     508        555    1095707                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform                     525        548    1109329                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform_mean                524        554    5492711                       // NOLINT(whitespace/line_length)
// BM_DoubleDiscreteCosineTransform_stddev               10          8    5492711                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    611        640     950628                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    582        632     963177                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    628        640     950133                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    583        632     962308                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform                    611        633     960732                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform_mean               603        635    4786978                       // NOLINT(whitespace/line_length)
// BM_TestTypeiscreteCosineTransform_stddev              18          4    4786978                       // NOLINT(whitespace/line_length)

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
