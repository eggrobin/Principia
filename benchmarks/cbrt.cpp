
#include "benchmarks/cbrt.hpp"

#include <string>

#include "glog/logging.h"
#include "quantities/quantities.hpp"

#include "mpirxx.h"

#if 0
#include "Intel/IACA 2.1/iacaMarks.h"
#define IACA_VOLATILE volatile
#else
#define IACA_VOLATILE
#define IACA_VC64_START
#define IACA_VC64_END
#endif

namespace principia {

using numerics::to_double;
using numerics::to_integer;

namespace numerics {

CubeRootRegistry& CubeRootRegistry::Instance() {
  static CubeRootRegistry* registry = new CubeRootRegistry();
  return *registry;
}

SingleParameterScalarFunction* CubeRootRegistry::Register(
    std::string const& name,
    SingleParameterScalarFunction* cbrt) {
  CHECK(methods_.emplace(name, cbrt).second) << name;
  return cbrt;
}

std::map<std::string, SingleParameterScalarFunction*> const&
CubeRootRegistry::methods() const {
  return methods_;
}

RoundedReal correct_cube_root(double const y) {
  RoundedReal result;
  int exponent;
  double y_mantissa = std::frexp(y, &exponent);
  while (exponent % 3 != 0) {
    y_mantissa *= 2;
    --exponent;
  }
  mpz_class const y_integer = y_mantissa * 0x1p300;
  mpz_class x_integer;
  mpz_root(x_integer.get_mpz_t(), y_integer.get_mpz_t(), 3);
  double x_towards_0 = x_integer.get_d();
  int x_inflated_exponent;
  std::frexp(x_towards_0, &x_inflated_exponent);
  double round_bit = std::ldexp(1, x_inflated_exponent - 54);
  mpz_class x_residual = x_integer - x_towards_0;
  if (x_residual == round_bit) {
    LOG(FATAL) << "truncation yields a tie for Y=" << to_integer(y);
  }
  double const towards_0 = x_towards_0 * 0x1p-100 * std::ldexp(1, exponent / 3);
  double const towards_0_ulps =
      mpz_class(x_towards_0 - x_integer).get_d() / (2 * round_bit);
  double const away_from_0 =
      (x_towards_0 + 2 * round_bit) * 0x1p-100 * std::ldexp(1, exponent / 3);
  double const away_from_0_ulps =
      mpz_class((x_towards_0 + 2 * round_bit) - x_integer).get_d() /
      (2 * round_bit);
  if (x_residual > round_bit) {
    result.nearest_rounding = away_from_0;
    result.furthest_rounding = towards_0;
    result.nearest_ulps = away_from_0_ulps;
    CHECK_LT(std::abs(result.nearest_ulps), 0.5);
  } else {
    result.nearest_rounding = towards_0;
    result.furthest_rounding = away_from_0;
    result.nearest_ulps = towards_0_ulps;
    CHECK_LT(std::abs(result.nearest_ulps), 0.5);
  }
  return result;
}

namespace egg {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
double cbrt(double const IACA_VOLATILE input) {
  IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  // Approximate ∛y with an error below 3,2 %.
  std::uint64_t const Y = to_integer(y);
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q³ = q * q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ = q - (q³ - y) * q / (2 * q³ + y);
  std::uint64_t const Ξ = to_integer(ξ) & 0xFFFF'FFF0'0000'0000;
  double const x = to_double(Ξ);
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
}  // namespace egg_halley

PRINCIPIA_REGISTER_CBRT(egg);

namespace egg_signed {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
double cbrt(double const IACA_VOLATILE input) {
  IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  // Approximate ∛y with an error below 3,2 %.
  std::uint64_t Y = to_integer(y);
  std::uint64_t const sign = 0x8000'0000'0000'0000 & Y;
  Y &= ~0x8000'0000'0000'0000;
  double const abs_y = to_double(Y);
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q³ = q * q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ = q - (q³ - abs_y) * q / (2 * q³ + abs_y);
  std::uint64_t const Ξ = to_integer(ξ) & 0xFFFF'FFF0'0000'0000;
  double const x = to_double(Ξ);
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const numerator =
      x * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x - numerator / denominator;
  double const IACA_VOLATILE signed_result =
      to_double(to_integer(result) | sign);
  IACA_VC64_END
  return signed_result;
}
}  // namespace egg_signed

namespace egg_signed_intrinsic {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128i sign_bit = _mm_cvtsi64_si128(0x8000'0000'0000'0000);
static const __m128i sixteen_bits_of_mantissa =
    _mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000);
double cbrt(double const IACA_VOLATILE input) {
  IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128i Y_0 = _mm_castpd_si128(_mm_set_sd(y));
  __m128i const sign = _mm_and_si128(sign_bit, Y_0);
  Y_0 = _mm_andnot_si128(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(_mm_castsi128_pd(Y_0));
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(Y_0);
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q³ = q * q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ = q - (q³ - abs_y) * q / (2 * q³ + abs_y);
  double const x = _mm_cvtsd_f64(_mm_castsi128_pd(_mm_and_si128(
      _mm_castpd_si128(_mm_set_sd(ξ)), sixteen_bits_of_mantissa)));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const numerator =
      x * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x - numerator / denominator;
  double const IACA_VOLATILE signed_result = _mm_cvtsd_f64(
      _mm_castsi128_pd(
      _mm_or_si128(
      _mm_castpd_si128(_mm_set_sd(result)), sign)));
  IACA_VC64_END
  return signed_result;
}
}  // namespace egg_signed_intrinsic

PRINCIPIA_REGISTER_CBRT(egg_signed_intrinsic);

#if PRINCIPIA_BENCHMARKS
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
                 quantities::DebugString(cbrt(2)) + u8"; -∛-2 = " +
                 quantities::DebugString(-cbrt(-2)));
}

void BM_EggCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &egg::cbrt);
}

void BM_EggSignedCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &egg_signed::cbrt);
}

void BM_EggSignedIntrinsicCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &egg_signed_intrinsic::cbrt);
}

BENCHMARK(BM_EggCbrt);
BENCHMARK(BM_EggSignedCbrt);
BENCHMARK(BM_EggSignedIntrinsicCbrt);
#endif

}  // namespace numerics
}  // namespace principia
