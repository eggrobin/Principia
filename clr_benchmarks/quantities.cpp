
// .\Release\clr_benchmarks.exe --benchmark_filter=Cosine --benchmark_min_time=3 --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/02-16:48:52
// Benchmark                                           Time(ns)    CPU(ns) Iterations                       // NOLINT(whitespace/line_length)
// ----------------------------------------------------------------------------------                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            546421     560225       1086                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            568150     538411       1130                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            551998     543218       1120                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            521166     542733       1121                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            551803     541767       1123                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform_mean       547950     545165       5580                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform_stddev      15280       7593       5580                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1833       1964     309854                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1855       1856     327738                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1851       1853     328387                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1767       1848     329243                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1903       1876     324343                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform_mean               1842       1878    1619565                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform_stddev               44         43    1619565                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          1848       1876     324280                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          1737       1845     329722                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          1846       1848     329205                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          1814       1867     325814                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          1850       1852     328600                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform_mean                     1819       1858    1637621                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform_stddev                     43         12    1637621                       // NOLINT(whitespace/line_length)

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
