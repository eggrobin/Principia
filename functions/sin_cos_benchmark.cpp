#include <algorithm>
#include <limits>
#include <random>

#include "boost/multiprecision/cpp_int.hpp"
#include "functions/benchmarking.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "numerics/next.hpp"
#include "numerics/sin_cos.hpp"
#include "quantities/numbers.hpp"
#include "testing_utilities/almost_equals.hpp"

// This test lives in `functions` to avoid pulling `boost` into `numerics`.
namespace principia {
namespace numerics {
namespace _sin_cos {

using namespace functions::_benchmarking;

class SinCosBenchmark : public ::testing::Test {
 protected:
  std::mt19937_64 random_{42};
  std::uniform_real_distribution<> uniformly_at_{-2 * π, 2 * π};
};

TEST_F(SinCosBenchmark, StdSinLatency) {
  // Note that we need to wrap std::sin in a lambda because of differences in
  // calling convention.
  std::cout << "std::sin latency: "
            << BenchmarkFunctionLatency(
                   &std::sin, [this]() { return uniformly_at_(random_); }, 10000)
                   .ToGUMString()
            << " cycles\n";
}

TEST_F(SinCosBenchmark, PrincipiaSinLatency) {
  std::cout << "Principia Sin latency: "
            << BenchmarkFunctionLatency(
                   &Sin, [this]() { return uniformly_at_(random_); }, 10000)
                   .ToGUMString()
            << " cycles\n";
}

TEST_F(SinCosBenchmark, StdSinThroughput) {
  std::cout << "std::sin reciprocal throughput: "
            << BenchmarkFunctionThroughput(
                   &std::sin, [this]() { return uniformly_at_(random_); }, 10000)
                   .ToGUMString()
            << " cycles\n";
}

TEST_F(SinCosBenchmark, PrincipiaSinThroughput) {
  std::cout << "Principia Sin reciprocal throughput: "
            << BenchmarkFunctionThroughput(
                   &Sin, [this]() { return uniformly_at_(random_); }, 10000)
                   .ToGUMString()
            << " cycles\n";
}

}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
