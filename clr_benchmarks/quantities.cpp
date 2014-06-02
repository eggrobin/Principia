
// .\Release\clr_benchmarks.exe --benchmark_filter=Cosine --benchmark_min_time=3 --benchmark_repetitions=5  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 2672 MHz CPU
// 2014/06/02-16:56:24
// Benchmark                                           Time(ns)    CPU(ns) Iterations                       // NOLINT(whitespace/line_length)
// ----------------------------------------------------------------------------------                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            536440     552091       1102                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            580556     562816       1081                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            540986     555620       1095                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            565101     555620       1095                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform            565374     556128       1094                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform_mean       557604     556433       5467                       // NOLINT(whitespace/line_length)
// BM_CLR_DimensionfulDiscreteCosineTransform_stddev      16526       3484       5467                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1781       1894     321211                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1834       1836     331399                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1834       1836     331404                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1755       1835     331470                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform                    1887       1863     326545                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform_mean               1818       1853    1642029                       // NOLINT(whitespace/line_length)
// BM_CLR_DoubleDiscreteCosineTransform_stddev               46         23    1642029                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          2677       2816     216052                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          2792       2778     218987                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          2790       2776     219152                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          2697       2802     217128                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform                          2788       2774     219328                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform_mean                     2749       2789    1090647                       // NOLINT(whitespace/line_length)
// BM_CLR_TestTypeCosineTransform_stddev                     51         17    1090647                       // NOLINT(whitespace/line_length)

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
