
#include "benchmarks/cbrt.hpp"

#include <string>

#include "glog/logging.h"
#include "quantities/quantities.hpp"
#include "numerics/fma.hpp"

#include "mpirxx.h"

#if 1
#include "Intel/IACA 2.1/iacaMarks.h"
#define IACA_FUNCTION_DOUBLE(arg) \
  (volatile double volatilizer) { \
    IACA_VC64_START;              \
    for (double iaca_result, arg = volatilizer;; arg = iaca_result)
#define IACA_RETURN(result)  \
  iaca_result = result;      \
  volatilizer = iaca_result; \
  IACA_VC64_END;             \
  }                          \
  return volatilizer
#else
#define IACA_FUNCTION_DOUBLE(arg) (double const arg)
#define IACA_RETURN(result) return result
#endif

namespace principia {

using numerics::to_double;
using numerics::to_integer;
using numerics::FusedMultiplyAdd;

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

namespace lagny_rational {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q³ = q * q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ = q - (q³ - abs_y) * q / (2 * q³ + abs_y);
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_rational

PRINCIPIA_REGISTER_CBRT(lagny_rational);

namespace lagny_rational_together {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q² = q * q;
  double const q³ = q² * q;
  double const q4 = q² * q²;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ = (q4 + 2 * y * q) / (2 * q³ + y);
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_rational_together

PRINCIPIA_REGISTER_CBRT(lagny_rational_together);

namespace lagny_rational_weighted {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q³ = q * q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ = q - (q³ - 1.01 * abs_y) * q / (2 * q³ + 1.02 * abs_y);
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_rational_weighted

PRINCIPIA_REGISTER_CBRT(lagny_rational_weighted);

namespace lagny_irrational_preinvert {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q³ = q * q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  constexpr double one_twelfth = 1 / 12.0;
  double const inverse = one_twelfth / q;
  __m128d const ρ_0 = _mm_set_sd((4 * abs_y - q³) * inverse);
  double const ξ = 0.5 * q + _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0));
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_irrational_preinvert

PRINCIPIA_REGISTER_CBRT(lagny_irrational_preinvert);

namespace lagny_irrational_expanded {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q² = q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  constexpr double one_twelfth = 1 / 12.0;
  __m128d const ρ_0 = _mm_set_sd(y / (3 * q) - q² * one_twelfth);
  double const ξ = 0.5 * q + _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0));
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_irrational_expanded

PRINCIPIA_REGISTER_CBRT(lagny_irrational_expanded);

namespace lagny_irrational_expanded_preinvert {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q² = q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  constexpr double one_twelfth = 1 / 12.0;
  constexpr double one_third = 1 / 3.0;
  double const inverse = one_third / q;
  __m128d const ρ_0 = _mm_set_sd(y * inverse - q² * one_twelfth);
  double const ξ = 0.5 * q + _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0));
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_irrational_expanded_preinvert

PRINCIPIA_REGISTER_CBRT(lagny_irrational_expanded_preinvert);

namespace lagny_irrational {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q³ = q * q * q;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  __m128d const ρ_0 = _mm_set_sd((4 * abs_y - q³) / (12 * q));
  double const ξ = 0.5 * q + _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0));
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_irrational

PRINCIPIA_REGISTER_CBRT(lagny_irrational);

namespace lagny_irrational_extracted_denominator {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  constexpr double sqrt_one_twelfth = 0.28867513459481288225457439025097872782380087563506;
  double const q² = q * q;
  double const q⁴ = q² * q²;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  __m128d const ρ_0 = _mm_set_sd(4 * abs_y * q  - q⁴);
  double const ξ = 0.5 * q + sqrt_one_twelfth / q * _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0));
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_irrational_extracted_denominator

PRINCIPIA_REGISTER_CBRT(lagny_irrational_extracted_denominator);

namespace lagny_canon_irrational_extracted_denominator {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q² = q * q;
  double const q⁴ = q² * q²;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  __m128d const ρ_0 = _mm_set_sd(4.00298737793169718250674332690180421 * abs_y * q  - q⁴);
  double const ξ = 0.499999938108574047751429172928306529 * q + 0.288531511562316719053845144194384063 / q * _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0));
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_canon_irrational_extracted_denominator

PRINCIPIA_REGISTER_CBRT(lagny_canon_irrational_extracted_denominator);

namespace lagny_canon_irrational_extracted_denominator5 {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt IACA_FUNCTION_DOUBLE(y) {
  __m128d sixteen_bits_of_mantissa{};
  sixteen_bits_of_mantissa = _mm_move_sd(
      sixteen_bits_of_mantissa,
      lagny_canon_irrational_extracted_denominator5::sixteen_bits_of_mantissa);
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q² = q * q;
  double const q⁴ = q² * q²;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  __m128d const ρ_0 = _mm_set_sd(4.00298737793169718250674332690180421 * abs_y * q  - q⁴);
  double const ξ = 0.499999938108574047751429172928306529 * q + 0.288531511562316719053845144194384063 / q * _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0));
  double const x = _mm_cvtsd_f64(_mm_and_pd(sixteen_bits_of_mantissa, _mm_set_sd(ξ)));
  // One round of 6th order Householder.
  double const x² = x * x;
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator = (x³ - abs_y) * ((10 * x³ + 16 * abs_y) * x³ + y²);
  double const denominator = x² * ((15 * x³ + 51 * abs_y) * x³ + 15 * y²);
  IACA_RETURN(x_sign_y - numerator / denominator);
}
}  // namespace lagny_canon_irrational_extracted_denominator5

PRINCIPIA_REGISTER_CBRT(lagny_canon_irrational_extracted_denominator5);

namespace lagny_canon_irrational_extracted_denominator_nearest {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q² = q * q;
  double const q⁴ = q² * q²;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  __m128d const ρ_0 = _mm_set_sd(4.00298737793169718250674332690180421 * abs_y * q  - q⁴);
  double const ξ = 0.499999938108574047751429172928306529 * q + 0.288531511562316719053845144194384063 / q * _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0));
  double const c = ξ * (0x1p36 + 1);
  double const x = (ξ - c) + c;
  // One round of 6th order Householder.
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const Δ = numerator / denominator;
  double const result = x_sign_y - Δ;
  double const residual = (x_sign_y - result) - Δ;
  double offset_from_halfway;
  if (residual > 0) {
    double const result_ulp_above =
        _mm_cvtsd_f64(_mm_castsi128_pd(_mm_add_epi64(
            _mm_castpd_si128(_mm_set_sd(result)), _mm_cvtsi64_si128(1)))) -
        result;
    offset_from_halfway = residual - 0.5 * result_ulp_above;
  } else {
    double const result_ulp_below =
        _mm_cvtsd_f64(_mm_castsi128_pd(_mm_sub_epi64(
            _mm_castpd_si128(_mm_set_sd(result)), _mm_cvtsi64_si128(1)))) -
        result;
    offset_from_halfway = residual - 0.5 * result_ulp_below;
  }
  double const distance_from_halfway =
      _mm_cvtsd_f64(_mm_andnot_pd(sign_bit, _mm_set_sd(offset_from_halfway)));
  possible_misrounding = distance_from_halfway < 0.000124 * 0x1p-53 * result;
  //IACA_VC64_END
  return result;
}
}  // namespace lagny_canon_irrational_extracted_denominator_nearest

PRINCIPIA_REGISTER_CBRT(lagny_canon_irrational_extracted_denominator_nearest);

namespace lagny_canon_irrational_extracted_denominator5_nearest {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d sixteen_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFF0'0000'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt IACA_FUNCTION_DOUBLE(y) {
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  double const q = to_double(Q);
  double const q² = q * q;
  double const q⁴ = q² * q²;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  __m128d const ρ_0 = _mm_set_sd(4.00298737793169718250674332690180421 * abs_y * q  - q⁴);
  double const ξ = 0.499999938108574047751429172928306529 * q + 0.288531511562316719053845144194384063 / q * _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0));
  double const c = ξ * (0x1p36 + 1);
  double const x = (ξ - c) + c;
  // One round of 6th order Householder.
  double const x² = x * x;
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator = (x³ - abs_y) * ((10 * x³ + 16 * abs_y) * x³ + y²);
  double const denominator = x² * ((15 * x³ + 51 * abs_y) * x³ + 15 * y²);
  double const Δ = numerator / denominator;
  IACA_RETURN(x_sign_y - Δ);
}
}  // namespace lagny_canon_irrational_extracted_denominator5_nearest

PRINCIPIA_REGISTER_CBRT(lagny_canon_irrational_extracted_denominator5_nearest);

namespace r5dr4_fma {
constexpr std::uint64_t C = 0x2A9F7893782DA1CE;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d twenty_five_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFFF'FA00'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
double cbrt(double const input) {
  //IACA_VC64_START
  double const y = input;
  // NOTE(egg): this needs rescaling and special handling of subnormal numbers.
  __m128d Y_0 = _mm_set_sd(y);
  __m128d const sign = _mm_and_pd(sign_bit, Y_0);
  Y_0 = _mm_andnot_pd(sign_bit, Y_0);
  double const abs_y = _mm_cvtsd_f64(Y_0);
  // Approximate ∛y with an error below 3,2 %.  I see no way of doing this with
  // SSE2 intrinsics, so we pay two cycles to move from the xmms to the r*xs and
  // back.
  std::uint64_t const Y = _mm_cvtsi128_si64(_mm_castpd_si128(Y_0));
  std::uint64_t const Q = C + Y / 3;
  // ⁰¹²³⁴⁵⁶⁷⁸⁹
  double const q = to_double(Q);
  double const y² = y * y;
  double const y³ = y² * y;
  double const q² = q * q;
  double const q³ = q² * q;
  double const q⁵ = q³ * q²;
  double const q⁶ = q³ * q³;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ =
      FusedMultiplyAdd(q⁶,
                       (FusedMultiplyAdd(5 * q, q², 45 * y)),
                       (FusedMultiplyAdd(30 * q, q², y)) * y²) /
      FusedMultiplyAdd(q⁵, FusedMultiplyAdd(15 * q, q², 51 * y), 15 * y² * q²);
  double const x =
      _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), twenty_five_bits_of_mantissa));
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const result = x_sign_y - numerator / denominator;
  //IACA_VC64_END
  return result;
}
}  // namespace r5dr4_fma

#if PRINCIPIA_BENCHMARKS
void BenchmarkCbrt(benchmark::State& state, double (*cbrt)(double)) {
  double total = 0;
  double total_cycles = 0;
  int iterations = 0;
  std::int64_t n = benchmark::CPUInfo::Get().cycles_per_second / 100;
  while (state.KeepRunning()) {
    double x = 1000;
    auto const start = __rdtsc();
    for (std::int64_t i = 0; i < n; ++i) {
      x = cbrt(x);
    }
    auto const stop = __rdtsc();
    total_cycles += stop - start;
    total += x;
    ++iterations;
  }
  state.SetLabel(std::to_string(total_cycles / (n * iterations)) + " cycles " +
                 quantities::DebugString(total / iterations, 3) + u8"; ∛2 = " +
                 quantities::DebugString(cbrt(2)) + u8"; ∛-2 = " +
                 quantities::DebugString(cbrt(-2)) + u8"; 2⁻³⁴⁰ ∛2¹⁰²¹ = " +
                 quantities::DebugString(0x1p-340 * cbrt(0x1p1021)) +
                 u8"; 2³⁴¹ ∛2⁻¹⁰²² = " +
                 quantities::DebugString(0x1p341 * cbrt(0x1p-1022)) +
                 u8"; 2³⁵⁸ ∛2⁻¹⁰⁷³ = " +
                 quantities::DebugString(0x1p358 * cbrt(0x1p-1073)));
  /*LOG(ERROR) << quantities::DebugString(cbrt(egg_scaling::big));
  LOG(ERROR) << quantities::DebugString(
      cbrt(egg_scaling::big * egg_scaling::big_σ⁻³) * egg_scaling::big_σ);
  LOG(ERROR) << quantities::DebugString(cbrt(egg_scaling::smol));
  LOG(ERROR) << quantities::DebugString(
      cbrt(egg_scaling::smol * egg_scaling::smol_σ⁻³) * egg_scaling::smol_σ);*/
}

void BM_NoCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, [](double x) { return x; });
}

void BM_StdCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, [](double x) { return std::cbrt(x); });
}

void BM_StdSin(benchmark::State& state) {
  BenchmarkCbrt(state, [](double x) { return std::sin(x)+std::cos(x); });
}

void BM_LagnyRationalCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_rational::cbrt);
}

void BM_LagnyRationalTogetherCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_rational_together::cbrt);
}

void BM_LagnyRationalWeightedCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_rational_weighted::cbrt);
}

void BM_LagnyIrrationalCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_irrational::cbrt);
}

void BM_LagnyIrrationalExtractedDenominatorCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_irrational_extracted_denominator::cbrt);
}

void BM_LagnyCanonIrrationalExtractedDenominatorCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_canon_irrational_extracted_denominator::cbrt);
}

void BM_LagnyCanonIrrationalExtractedDenominator5Cbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_canon_irrational_extracted_denominator5::cbrt);
}

void BM_LagnyCanonIrrationalExtractedDenominatorNearestCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_canon_irrational_extracted_denominator_nearest::cbrt);
}

void BM_LagnyCanonIrrationalExtractedDenominator5NearestCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_canon_irrational_extracted_denominator5_nearest::cbrt);
}

void BM_LagnyIrrationalPreinvertCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_irrational_preinvert::cbrt);
}

void BM_LagnyIrrationalExpandedCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_irrational_expanded::cbrt);
}

void BM_LagnyIrrationalExpandedPreinvertCbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &lagny_irrational_expanded_preinvert::cbrt);
}

void BM_R5DR4FMACbrt(benchmark::State& state) {
  BenchmarkCbrt(state, &r5dr4_fma::cbrt);
}

BENCHMARK(BM_NoCbrt);
BENCHMARK(BM_LagnyRationalCbrt);
BENCHMARK(BM_LagnyRationalTogetherCbrt);
BENCHMARK(BM_LagnyRationalWeightedCbrt);
BENCHMARK(BM_LagnyIrrationalCbrt);
BENCHMARK(BM_LagnyIrrationalExpandedCbrt);
BENCHMARK(BM_LagnyIrrationalPreinvertCbrt);
BENCHMARK(BM_LagnyIrrationalExpandedPreinvertCbrt);
BENCHMARK(BM_LagnyIrrationalExtractedDenominatorCbrt);
BENCHMARK(BM_LagnyCanonIrrationalExtractedDenominatorCbrt);
BENCHMARK(BM_LagnyCanonIrrationalExtractedDenominator5Cbrt);
BENCHMARK(BM_R5DR4FMACbrt);
BENCHMARK(BM_LagnyCanonIrrationalExtractedDenominatorNearestCbrt);
BENCHMARK(BM_LagnyCanonIrrationalExtractedDenominator5NearestCbrt);
BENCHMARK(BM_StdCbrt);
BENCHMARK(BM_StdSin);
#endif

}  // namespace numerics
}  // namespace principia
