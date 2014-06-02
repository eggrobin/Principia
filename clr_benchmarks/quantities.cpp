
// .\Release\clr_benchmarks.exe --benchmark_filter=Cosine --benchmark_min_time=3 --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/02-13:41:56
// Benchmark                                           Time(ns)    CPU(ns) Iterations                       // NOLINT(whitespace/line_length)
// ----------------------------------------------------------------------------------                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            554878     552592       1101                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            566269     549101       1108                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            520682     543218       1120                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            551431     542733       1121                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            536841     542733       1121                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform_mean       545946     546046       5571                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform_stddev      15779       4043       5571                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                      76        107    5660370                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                      73        104    5827506                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                      67        105    5800195                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                      73        105    5818900                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                      71        106    5747258                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform_mean                 72        105   28854229                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform_stddev                3          1   28854229                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                            73        106    5761202                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                            72        104    5854591                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                            69        105    5771305                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                            73        104    5843699                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                            72        104    5862314                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform_mean                       72        105   29093111                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform_stddev                      1          1   29093111                       // NOLINT(whitespace/line_length)

#include "benchmark/benchmark.h"

using principia::clr_benchmarks_adapter::QuantitiesCLRBenchmark;

namespace principia {
namespace clr_benchmarks {

static void BM_CLR_DimensionfulDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    QuantitiesCLRBenchmark::DimensionfulDiscreteCosineTransform();
  }
}

BENCHMARK(BM_CLR_DimensionfulDiscreteCosineTransform);

static void BM_CLR_DoubleDiscreteCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    QuantitiesCLRBenchmark::DoubleDiscreteCosineTransform();
  }
}

BENCHMARK(BM_CLR_DoubleDiscreteCosineTransform);

static void BM_CLR_TestTypeCosineTransform(
    benchmark::State& state) {  // NOLINT(runtime/references)
  while (state.KeepRunning()) {
    QuantitiesCLRBenchmark::TestTypeDiscreteCosineTransform();
  }
}

BENCHMARK(BM_CLR_TestTypeCosineTransform);

}  // namespace clr_benchmarks
}  // namespace principia
