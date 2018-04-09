﻿
#include <string>

#include "glog/logging.h"
#if !NO_BENCHMARK
#include "benchmark/benchmark.h"
#endif
#include "quantities/quantities.hpp"

#include "mpirxx.h"

#include <emmintrin.h>
#if 0
#include "Intel/IACA 2.1/iacaMarks.h"
#define IACA_VOLATILE volatile
#else
#define IACA_VOLATILE
#define IACA_VC64_START
#define IACA_VC64_END
#endif

namespace principia {

namespace {
std::uint64_t to_integer(double x) {
  return _mm_cvtsi128_si64(_mm_castpd_si128(_mm_set_sd(x)));
}
double to_double(std::uint64_t x) {
  return _mm_cvtsd_f64(_mm_castsi128_pd(_mm_cvtsi64_si128(x)));
}
}  // namespace

namespace slow_correct {
double cbrt(double const y) {
  int exponent;
  double y_mantissa = std::frexp(y, &exponent);
  while (exponent % 3 != 0) {
    y_mantissa *= 2;
    --exponent;
  }
  mpz_class const y_integer = y_mantissa * 0x1p201;
  mpz_class x_integer;
  mpz_root(x_integer.get_mpz_t(), y_integer.get_mpz_t(), 3);
  double x_towards_0 = x_integer.get_d();
  int x_inflated_exponent;
  std::frexp(x_towards_0, &x_inflated_exponent);
  double round_bit = std::ldexp(1, x_inflated_exponent - 54);
  mpz_class x_residual = x_integer - x_towards_0;
  if (x_residual == round_bit) {
    LOG(FATAL) << "truncation yields a tie";
  }
  if (x_residual > round_bit) {
    return (x_towards_0 + 2 * round_bit) * 0x1p-67 *
           std::ldexp(1, -exponent / 3);
  } else {
    return x_towards_0 * 0x1p-67 * std::ldexp(1, -exponent / 3);
  }
}
}  // namespace slow_correct

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

namespace one_halley_iterate {
constexpr std::uint64_t C = 0x2a9f76253119d328;
double cbrt(double const y) {
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  std::uint64_t Y = to_integer(y);
  std::uint64_t Q = C + Y / 3;
  double x = to_double(Q);
  // One M = 3 iterate.
  double x³ = x * x * x;
  x = x - (x³ - y) * x / (2 * x³ + y);
  return x;
}
}  // namespace one_halley_iterate

namespace householder_order_10 {
constexpr std::uint64_t C = 0x2a9f76253119d328;
double cbrt(double const y) {
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  std::uint64_t Y = to_integer(y);
  std::uint64_t Q = C + Y / 3;
  double const x = to_double(Q);
  double const x³ = x * x * x;
  double const y² = y * y;
  double const y³ = y * y²;
  double const y⁴ = y² * y²;
  double const y⁵ = y² * y³;
  double const y⁶ = y³ * y³;
  double const numerator_big_factor =
      (y⁵ + x³ * (46 * y⁴ +
                  x³ * (256 * y³ + x³ * (323 * y² + x³ * (5 * x³ + 98 * y)))));
  double const numerator = 9 * x * (x³ - y) * numerator_big_factor;
  double const denominator =
      y⁶ + x³ * (210 * y⁵ +
                 x³ * (2850 * y⁴ +
                       x³ * (8350 * y³ +
                             x³ * (6765 * y² + x³* (55 * x³ + 1452 * y)))));
  return x - numerator / denominator;
}
}  // namespace householder_order_10

namespace householder_order_10_estrin {
constexpr std::uint64_t C = 0x2a9f76253119d328;
double cbrt(double const y) {
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  std::uint64_t Y = to_integer(y);
  std::uint64_t Q = C + Y / 3;
  double const x = to_double(Q & 0xFFFF'FFF0'0000'0000);
  double const x³ = x * x * x;
  double const y² = y * y;
  double const y⁴ = y² * y²;
  double const x⁶ = x³ * x³;
  double const x⁹ = x⁶ * x³;
  double const numerator_big_factor =
      ((5 * x³ + 98 * y) * x⁶ + (323 * x³ + 256 * y) * y²) * x⁶ +
      (46 * x³ + y) * y⁴;
  double const numerator = 9 * x * (x³ - y) * numerator_big_factor;
  double const denominator =
      ((55 * x³ + 1452 * y) * x⁶ + (6765 * x³ + 8350 * y) * y²) * x⁹ +
      ((2850 * x³ + 210 * y) * x³ + y²) * y⁴;
  return x - numerator / denominator;
}
}  // namespace householder_order_10_estrin

namespace kahan {
constexpr std::uint64_t C = 0x2a9f76253119d328;
double cbrt(double const y) {
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  std::uint64_t const Y = to_integer(y);
  std::uint64_t const Q = C + Y / 3;
  double x = to_double(Q);
  // One M = 3 iterate.
  double x³ = x * x * x;
  x = x - (x³ - y) * x / (2 * x³ + y);
  std::uint64_t const X = to_integer(x) & 0xFFFF'FFF0'0000'0000;
  x = to_double(X);
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
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
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

namespace egg {
constexpr std::uint64_t G = 0x553ef0ff289dd796;
double cbrt(double const IACA_VOLATILE input) {
  IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  // Approximate 1/∛y with an error below 3,5 %.
  std::uint64_t const Y = to_integer(y);
  std::uint64_t const R = G - Y / 3;
  double const r = to_double(R);
  // z = z₁z₂ is the approximation of 1/∛y by two rounds of Newton on r.
  // TODO(egg): error here.
  double const r³y = (r * r) * (r * y);
  double const r⁶y² = r³y * r³y;
  double const z₁ = -1.0 / 243.0 * r * (r³y - 4);
  double const z₂ = ((108 - 64 * r³y) + (48 - 12 * r³y + r⁶y²) * r⁶y²);
  // An approximation of ∛y [TODO(egg): error here].
  double const yz² = y * (z₁ * z₁) * (z₂ * z₂);
  std::uint64_t const X = to_integer(yz²) & 0xFFFF'FFF0'0000'0000;
  double const x = to_double(X);
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const numerator = x * (x³ - y) * ((5 * x³ + 17 * y) * x³ + 5 * y²);
  double const denominator = (7 * x³ + 42 * y) * x⁶ + (30 * x³ + 2 * y) * y²;
  double const IACA_VOLATILE result = x - numerator / denominator;
  IACA_VC64_END
  return result;
}
}  // namespace egg

#if !NO_BENCHMARK
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

void BM_OneHalleyIterateCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &one_halley_iterate::cbrt);
}

void BM_HouseholderOrder10EstrinCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &householder_order_10_estrin::cbrt);
}

void BM_HouseholderOrder10Cbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &householder_order_10::cbrt);
}

void BM_KahanCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &kahan::cbrt);
}

void BM_KahanNoDivCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &kahan_no_div::cbrt);
}

void BM_EggCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &egg::cbrt);
}

void BM_MicrosoftCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &std::cbrt);
}

void BM_SlowCorrectCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &slow_correct::cbrt);
}

BENCHMARK(BM_AtlasCbrt);
BENCHMARK(BM_OneHalleyIterateCbrt);
BENCHMARK(BM_HouseholderOrder10EstrinCbrt);
BENCHMARK(BM_HouseholderOrder10Cbrt);
BENCHMARK(BM_KahanCbrt);
BENCHMARK(BM_KahanNoDivCbrt);
BENCHMARK(BM_EggCbrt);
BENCHMARK(BM_MicrosoftCbrt);
BENCHMARK(BM_SlowCorrectCbrt);
#endif

}  // namespace principia
