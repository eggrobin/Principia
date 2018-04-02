
#include <string>

#include "glog/logging.h"
#include "benchmark/benchmark.h"
#include "quantities/quantities.hpp"

namespace principia {

namespace {
std::uint64_t to_integer(double x) {
  std::uint64_t result;
  std::memcpy(&result, &x, sizeof(x));
  return result;
}
double to_double(std::uint64_t x) {
  double result;
  std::memcpy(&result, &x, sizeof(x));
  return result;
}
}  // namespace

namespace atlas {
//  curl https://raw.githubusercontent.com/simonbyrne/apple-libm/4853bcad08357b0d7991d46bf95d578a909be27b/Source/ARM/cbrt.c `
//  | select -expandproperty content > `
// .\benchmarks\atlas_cbrt

// we ignore subnormal numbers in the other implementations, so let's make the
// comparison fair.
#define DAZ 0
#define AVOID_UINT64 0
#include "benchmarks/atlas_cbrt"
}  // namespace atlas

namespace kahan {
constexpr std::uint64_t C = 0x2a9f76253119d328;
double cbrt(double const y) {
  // NOTE(eggrobin): this needs rescaling and special handling of subnormal
  // numbers.
  std::uint64_t Y = to_integer(y);
  std::uint64_t Q = C + Y / 3;
  double x = to_double(Q);
  // One M = 3 iterate.
  double x³ = x * x * x;
  x = x - (x³ - y) * x / (2 * x³ + y);
  // One M = 4 iterate, Γ = -35/3.
  x³ = x * x * x;
  double const s = (x³ - y) / y;
  x = x - x * ((14.0 / 81.0 * s - 2.0 / 9.0) * s + 1.0 / 3.0) * s;
  return x;
}
}  // namespace kahan

namespace kahan_no_div {
constexpr std::uint64_t G  = 0x553ef0ff289dd796;
double cbrt(double const y) {
  // NOTE(eggrobin): this needs rescaling and special handling of subnormal
  // numbers.
  std::uint64_t Y = to_integer(y);
  std::uint64_t R = G - Y / 3;
  double z = to_double(R);
  double z²;
  double z⁴;
  z² = z * z;
  z⁴ = z² * z²;
  z = z + (1.0 / 3.0) * (z - z⁴ * y);
  z² = z * z;
  z⁴ = z² * z²;
  z = z + (1.0 / 3.0) * (z - z⁴ * y);
  z² = z * z;
  z⁴ = z² * z²;
  z = z + (1.0 / 3.0) * (z - z⁴ * y);
  z² = z * z;
  z⁴ = z² * z²;
  z = z + (1.0 / 3.0) * (z - z⁴ * y);
  double const x = z * z * y;
  return x;
}
}  // namespace kahan_no_div

void BenchmarkCbrt(benchmark::State& state, double (*cbrt)(double)) {
  double total = 0;
  int iterations = 0;
  while (state.KeepRunning()) {
    double x = 1000;
    for (int i = 0; i < 1000; ++i) {
      x = cbrt(x);
    }
    total += x;
    ++iterations;
  }
  state.SetLabel(quantities::DebugString(total / iterations) + u8"; ∛2 = " +
                 quantities::DebugString(cbrt(2)));
}

void BM_AtlasCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &atlas::cbrt);
}

void BM_KahanCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &kahan::cbrt);
}

void BM_KahanNoDivCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &kahan_no_div::cbrt);
}

void BM_MicrosoftCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &std::cbrt);
}

BENCHMARK(BM_AtlasCbrt);
BENCHMARK(BM_KahanCbrt);
BENCHMARK(BM_KahanNoDivCbrt);
BENCHMARK(BM_MicrosoftCbrt);

}  // namespace principia
