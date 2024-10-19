#pragma once

#include <functional>
#include <string>

namespace principia {
namespace functions {
namespace _benchmarking {

struct MeasurementResult {
  double value{};
  double standard_uncertainty{};
  std::string ToGUMString() const;
};

MeasurementResult BenchmarkFunctionThroughput(double (__cdecl *f)(double),
                                              std::function<double()> get_input,
                                              std::int64_t const samples);

MeasurementResult BenchmarkFunctionLatency(double (__cdecl *f)(double),
                                           std::function<double()> get_input,
                                           std::int64_t const samples);

}  // namespace _benchmarking
}  // namespace functions
}  // namespace principia