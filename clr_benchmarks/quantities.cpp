
// .\Release\clr_benchmarks.exe --benchmark_filter=Cosine --benchmark_min_time=3 --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/02-17:08:14
// Benchmark                                           Time(ns)    CPU(ns) Iterations                       // NOLINT(whitespace/line_length)
// ----------------------------------------------------------------------------------                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            545258     560225       1086                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            551852     542733       1121                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            551161     542249       1122                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            527524     541767       1123                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            567624     551590       1103                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform_mean       548637     547618       5555                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform_stddev      12935       7195       5555                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1843       1868     325719                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1733       1836     331308                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1837       1839     330779                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1807       1863     326568                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1836       1835     331584                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform_mean               1811       1848    1645958                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform_stddev               41         14    1645958                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                        605775     602380       1010                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                        576961     592986       1026                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                        619484     600596       1013                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                        577000     593565       1025                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                        603621     592409       1027                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform_mean                   596486     596357       5101                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform_stddev                  16892       4212       5101                       // NOLINT(whitespace/line_length)

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
