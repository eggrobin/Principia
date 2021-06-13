
#include "benchmarks/cbrt.hpp"

#include <array>
#include <random>
#include <string>

#include "glog/logging.h"
#include "quantities/quantities.hpp"
#include "quantities/elementary_functions.hpp"
#include "numerics/fma.hpp"
#include "numerics/double_precision.hpp"

#include "mpirxx.h"

#define NOIACA_FUNCTION_DOUBLE(arg) (double const arg)
#define NOIACA_RETURN(result) return result
#if 0
#include "Intel/IACA 2.1/iacaMarks.h"
#include <random>
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

std::array<double, 4> NievergeltQuadruplyCompensatedStep(
    DoublePrecision<double> b,
    DoublePrecision<double> d) {
  auto const& [b₁, b₂] = b;
  auto const& [d₁, d₂] = d;
  if (std::abs(d₂) >= std::abs(b₁)) {
    std::swap(b, d);
  }
  if (std::abs(b₂) >= std::abs(d₁)) {
    double const g₄ = d₂;
    auto const [w, g₃] = TwoSum(b₂, d₁);
    auto const [g₁, g₂] = TwoSum(b₁, w);
    return {g₁, g₂, g₃, g₄};
  } else {
    auto const [w, g₄] = TwoSum(b₂, d₂);
    auto const [u, v] = TwoSum(d₁, b₁);
    auto const [x, g₃] = TwoSum(v, w);
    auto const [g₁, g₂] = TwoSum(u, x);
    return {g₁, g₂, g₃, g₄};
  }
}

template<std::size_t n>
std::array<double, n> PriestNievergeltNormalize(std::array<double, n> const f) {
  std::array<double, n> s{};
  int k = 0;
  s[0] = f[0];
  for (int j = 1; j < n; ++j) {
    auto const [c, d] = TwoSum(s[k], f[j]);
    s[k] = c;
    if (d != 0) {
      int l = k - 1;
      k = k + 1;
      while (l >= 0) {
        auto const [cʹ, dʹ] = TwoSum(s[l], s[l + 1]);
        s[l] = cʹ;
        if (dʹ == 0) {
          k = k - 1;
        } else {
          s[l + 1] = dʹ;
        }
        l = l - 1;
      }
      s[k] = d;
    }
  }
  return s;
}

void ConsiderCorrection(double const r₀, double const r₁, double const τ) {
  double const r̃ = r₀ + 2 * r₁;
  possible_misrounding = r̃ != r₀ && std::abs(0.5 * (r̃ - r₀) - r₁) <= τ * r₀;
}

double CorrectLastBit(double const y, double const r₀, double const r₁, double const τ) {
  double const r̃ = r₀ + 2 * r₁;
  if (r̃ == r₀ || std::abs(0.5 * (r̃ - r₀) - r₁) > τ * r₀) {
    return r₀;
  }
  // TODO(egg): Handle negative y.
  CHECK_GT(y, 0);
  double const a = std::min(r₀, r̃);
  double const b = 0.5 * (std::max(r₀, r̃) - a);
  double const b² = b * b;
  DoublePrecision<double> const a² = TwoProduct(a, a);
  auto const& [a²₀, a²₁] = a²;
  DoublePrecision<double> const a³₀ = TwoProduct(a²₀, a);
  DoublePrecision<double> minus_a³₁ = TwoProduct(a²₁, -a);
  auto const& [a³₀₀, a³₀₁] = a³₀;
  // ρ₅₃ = y - a³ = y - a³₀ - a³₁ = y - a³₀₀ - a³₀₁ - a³₁;
  double const ρ₀ = y - a³₀₀;  // Exact.
  // ρ₅₃ = ρ₀ - a³₀₁ - a³₁;
  std::array<double, 4> const ρ₅₃ = PriestNievergeltNormalize(
      NievergeltQuadruplyCompensatedStep(TwoDifference(ρ₀, a³₀₁), minus_a³₁));
  CHECK_EQ(ρ₅₃[3], 0);
  std::array<double, 3> ρ₅₄{ρ₅₃[0], ρ₅₃[1], ρ₅₃[2]};
  for (double rhs :
       {2 * a²₀ * b, a²₀ * b, 2 * a²₁ * b, a²₁ * b, 2 * a * b², a * b²}) {
    auto const ρ = PriestNievergeltNormalize(NievergeltQuadruplyCompensatedStep(
        TwoSum(ρ₅₄[0], ρ₅₄[1]), TwoDifference(ρ₅₄[2], rhs)));
    CHECK_EQ(ρ[3], 0);
    ρ₅₄ = {ρ[0], ρ[1], ρ[2]};
  }
  bool const ρ₅₄_positive = ρ₅₄[0] > 0 || (ρ₅₄[0] == 0 && ρ₅₄[1] > 0) ||
                            (ρ₅₄[0] == 0 && ρ₅₄[1] == 0 && ρ₅₄[2] >= 0);
  return ρ₅₄_positive ? std::max(r₀, r̃) : a;
}

namespace fast_correct {
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
__declspec(noinline) double cbrt NOIACA_FUNCTION_DOUBLE(y) {
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
  double const ξ = (q4 + 2 * abs_y * q) / (2 * q³ + abs_y);
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
  double const Δ = numerator / denominator;
  double const r₀ = x_sign_y - Δ;
  double const r₁ = x_sign_y - r₀ - Δ;
  return CorrectLastBit(y, r₀, r₁, /*τ=*/std::numeric_limits<double>::infinity());
}
}  // namespace fast_correct

PRINCIPIA_REGISTER_CBRT(fast_correct);

namespace kahan {
constexpr std::uint64_t C = 0x2A9F76253119D328;
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
__declspec(noinline) double cbrt NOIACA_FUNCTION_DOUBLE(y) {
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
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ = q - ((q³ - y) * q) / (2 * q³ + y);
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));

  double const x³ = x * x * x;
  //double const s = (x³ - abs_y) / abs_y;
  //return x - x * ((14.0 / 81 * s - 2.0 / 9) * s + 1.0 / 3) * s;
  double const λ = 0x1.16B28F55D72D4p-2;
  double const z = (x³ - abs_y);
  double const λz = λ * z;
  return x -
         abs_y * z * x / (3 * (abs_y - λz) * (abs_y + λz) + 2 * (abs_y * z));
}
}  // namespace kahan

PRINCIPIA_REGISTER_CBRT(kahan);

namespace r3dr6 {
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
__declspec(noinline) double cbrt NOIACA_FUNCTION_DOUBLE(y) {
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
  double const ξ = (q4 + 2 * abs_y * q) / (2 * q³ + abs_y);
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
  double const Δ = numerator / denominator;
  double const r₀ = x_sign_y - Δ;
  return r₀;
}
}  // namespace r3dr6

PRINCIPIA_REGISTER_CBRT(r3dr6);

namespace r3dr5 {
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
__declspec(noinline) double cbrt NOIACA_FUNCTION_DOUBLE(y) {
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
  double const ξ = (q4 + 2 * abs_y * q) / (2 * q³ + abs_y);
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x² = x * x;
  double const x³ = x * x * x;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const x²_sign_y = x_sign_y * x;
  double const numerator = (x³ - abs_y) * ((10 * x³ + 16 * abs_y) * x³ + y²);
  double const denominator =
      x²_sign_y * ((15 * x³ + 51 * abs_y) * x³ + 15 * y²);
  NOIACA_RETURN(x_sign_y - numerator / denominator);
}
}  // namespace r3dr5

PRINCIPIA_REGISTER_CBRT(r3dr5);

namespace i3tdr5 {
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
__declspec(noinline) double cbrt NOIACA_FUNCTION_DOUBLE(y) {
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
  __m128d const ρ_0 = _mm_set_sd(0x1.0030F1F8A11DAp2 * abs_y * q - q⁴);
  double inverse = 0x1.2774CDF81A35Ep-2 / q;
  double const ξ =
      (0x1.BBA02BAFEA9B7p0 * q² + _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0))) * inverse;
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  // One round of 6th order Householder.
  double const x² = x * x;
  double const x³ = x * x * x;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const x²_sign_y = x_sign_y * x;
  double const numerator = (x³ - abs_y) * ((10 * x³ + 16 * abs_y) * x³ + y²);
  double const denominator = x²_sign_y * ((15 * x³ + 51 * abs_y) * x³ + 15 * y²);
  double const Δ = numerator / denominator;
  double const r₀ = x_sign_y - Δ;
  double const r₁ = x_sign_y - r₀ - Δ;
  ConsiderCorrection(r₀, r₁, /*τ=*/0x1.7C73DBBD9FA60p-66);
  return r₀;
}
}  // namespace i3tdr5

PRINCIPIA_REGISTER_CBRT(i3tdr5);

namespace i3tdr6 {
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
__declspec(noinline) double cbrt NOIACA_FUNCTION_DOUBLE(y) {
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
  __m128d const ρ_0 = _mm_set_sd(0x1.0030F1F8A11DAp2 * abs_y * q - q⁴);
  double inverse = 0x1.2774CDF81A35Ep-2 / q;
  double const ξ =
      (0x1.BBA02BAFEA9B7p0 * q² + _mm_cvtsd_f64(_mm_sqrt_sd(ρ_0, ρ_0))) * inverse;
  double const x = _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), sixteen_bits_of_mantissa));
  double const x³ = x * x * x;
  double const x⁶ = x³ * x³;
  double const y² = y * y;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator =
      x_sign_y * (x³ - abs_y) * ((5 * x³ + 17 * abs_y) * x³ + 5 * y²);
  double const denominator =
      (7 * x³ + 42 * abs_y) * x⁶ + (30 * x³ + 2 * abs_y) * y²;
  double const Δ = numerator / denominator;
  double const r₀ = x_sign_y - Δ;
  double const r₁ = x_sign_y - r₀ - Δ;
  ConsiderCorrection(r₀, r₁, /*τ=*/0x1.AC20CF34393E1p-66);
  return r₀;
}
}  // namespace i3tdr6

PRINCIPIA_REGISTER_CBRT(i3tdr6);
namespace r5dr4_fma {
constexpr std::uint64_t C = 0x2A9F76253119D328;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d twenty_five_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFFF'F800'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
__declspec(noinline) double cbrt IACA_FUNCTION_DOUBLE(y) {
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
  double const y³ = y² * abs_y;
  double const q² = q * q;
  double const q³ = q² * q;
  double const q⁵ = q³ * q²;
  double const q⁶ = q³ * q³;
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ =
      FusedMultiplyAdd(q⁶,
                       (FusedMultiplyAdd(5 * q, q², 45 * abs_y)),
                       (FusedMultiplyAdd(30 * q, q², abs_y)) * y²) /
      FusedMultiplyAdd(
          q⁵, FusedMultiplyAdd(15 * q, q², 51 * abs_y), 15 * y² * q²);
  double const x =
      _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), twenty_five_bits_of_mantissa));
  double const x² = x * x;
  DCHECK_EQ(FusedMultiplySubtract(x, x, x²), 0);
  double const x³ = x² * x;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator = x_sign_y * FusedMultiplySubtract(x², x, abs_y);
  double const denominator =
      FusedMultiplyAdd(x³, FusedMultiplyAdd(10 * x, x², 16 * abs_y), y²);
  double const result =
      FusedNegatedMultiplyAdd(FusedMultiplyAdd(6 * x, x², 3 * abs_y),
                              numerator / denominator,
                              x_sign_y);
  IACA_RETURN(result);
}
}  // namespace r5dr4_fma

PRINCIPIA_REGISTER_CBRT(r5dr4_fma);

namespace i5dr4_fma {
constexpr std::uint64_t C = 0x2A9F76253119D328;
static const __m128d sign_bit =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
static const __m128d twenty_five_bits_of_mantissa =
    _mm_castsi128_pd(_mm_cvtsi64_si128(0xFFFF'FFFF'F800'0000));
// NOTE(egg): the σs do not rescale enough to put the least normal or greatest
// finite magnitudes inside the non-rescaling range; for very small and very
// large values, rescaling occurs twice.
constexpr double smol = 0x1p-225;
constexpr double smol_σ = 0x1p-154;
constexpr double smol_σ⁻³ = 1 / (smol_σ * smol_σ * smol_σ);
constexpr double big = 0x1p237;
constexpr double big_σ = 0x1p154;
constexpr double big_σ⁻³ = 1 / (big_σ * big_σ * big_σ);
__declspec(noinline) double cbrt IACA_FUNCTION_DOUBLE(y) {
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
  double const inverse =
      q / FusedMultiplySubtract(
              0x1.4A7E9CB8A3491p2 * q, q², 0x1.08654A2D4F6DBp-1 * abs_y);
  // An approximation of ∛y with a relative error below 2⁻¹⁵.
  double const ξ = FusedMultiplyAdd(
      inverse,
      quantities::Sqrt(FusedMultiplySubtract(
          q³, FusedNegatedMultiplyAdd(q², q, 118.0 / 5 * abs_y), y²)),
      FusedMultiplySubtract(q², q, abs_y) * 0x1.4A7E9CB8A3491p0 * inverse);
  double const x =
      _mm_cvtsd_f64(_mm_and_pd(_mm_set_sd(ξ), twenty_five_bits_of_mantissa));
  double const x² = x * x;
  DCHECK_EQ(FusedMultiplySubtract(x, x, x²), 0);
  double const x³ = x² * x;
  double const x_sign_y = _mm_cvtsd_f64(_mm_or_pd(_mm_set_sd(x), sign));
  double const numerator = x_sign_y * FusedMultiplySubtract(x², x, abs_y);
  double const denominator =
      FusedMultiplyAdd(x³, FusedMultiplyAdd(10 * x, x², 16 * abs_y), y²);
  double const Δ₁ = FusedMultiplyAdd(6 * x, x², 3 * abs_y);
  double const Δ₂ = numerator / denominator;
  double const r₀ = FusedNegatedMultiplyAdd(Δ₁, Δ₂, x_sign_y);
  double const r₁ = FusedNegatedMultiplyAdd(Δ₁, Δ₂, x_sign_y - r₀);
  // TODO(egg): ConsiderCorrection.
  return r₀;
}
}  // namespace i5dr4_fma

PRINCIPIA_REGISTER_CBRT(i5dr4_fma);

#if PRINCIPIA_BENCHMARKS

__declspec(noinline) void BenchmarkCbrtLatency(benchmark::State& state, double (*cbrt)(double)) {
  static double zero_cycles = quantities::Infinity<double>;
  double total = 0;
  double min_cycles = quantities::Infinity<double>;
  int iterations = 0;
  std::int64_t n = 1 << 16;
  constexpr std::uint64_t low = 0x3FF0000000000000;   // 1.
  constexpr std::uint64_t high = 0x4020000000000000;  // 8.
  std::linear_congruential_engine<std::uint64_t,
                                  6364136223846793005,
                                  1442695040888963407,
                                  0>
      rng(1729);
  while (state.KeepRunning()) {
    double x = 1000;
    auto const start = __rdtsc();
    for (std::int64_t i = 0; i < n; ++i) {
      rng.seed(i + principia::numerics::to_integer(x));
      std::uint64_t const Y = rng() % (high - low) + low;
      double const y = principia::numerics::to_double(Y);
      x = cbrt(y);
    }
    auto const stop = __rdtsc();
    min_cycles = std::min(min_cycles, static_cast<double>(stop - start));
    total += x;
    ++iterations;
  }
  if (cbrt(5) == 5) {
    zero_cycles = std::min(zero_cycles, min_cycles / n);
  }
  state.SetLabel(std::to_string(min_cycles / n - zero_cycles) +
                 " cycles " + quantities::DebugString(total, 3) + u8"; ∛-2 = " +
                 quantities::DebugString(cbrt(-2)));
}

__declspec(noinline) void BenchmarkCbrtThroughput(benchmark::State& state, double (*cbrt)(double)) {
  static double zero_cycles = quantities::Infinity<double>;
  double total = 0;
  double min_cycles = quantities::Infinity<double>;
  int iterations = 0;
  constexpr std::int64_t n = 64;
  constexpr std::uint64_t low = 0x3FF0000000000000;   // 1.
  constexpr std::uint64_t high = 0x4020000000000000;  // 8.
  std::linear_congruential_engine<std::uint64_t,
                                  6364136223846793005,
                                  1442695040888963407,
                                  0>
      rng(1729);
  std::array<double, static_cast<std::size_t>(n)> inputs;
  for (std::int64_t i = 0; i < n; ++i) {
    std::uint64_t const Y = rng() % (high - low) + low;
    double const y = principia::numerics::to_double(Y);
    inputs[i] = y;
  }
  while (state.KeepRunning()) {
    auto const start = __rdtsc();
    for (std::int64_t i = 0; i < n; ++i) {
      inputs[i] = cbrt(inputs[i]) * π;
    }
    auto const stop = __rdtsc();
    min_cycles = std::min(min_cycles, static_cast<double>(stop - start));
    ++iterations;
  }
  for (std::int64_t i = 0; i < n; ++i) {
    total += inputs[i] / n;
  }
  if (cbrt(5) == 5) {
    zero_cycles = std::min(zero_cycles, min_cycles / n);
  }
  state.SetLabel(std::to_string(min_cycles / n - zero_cycles) +
                 " cycles " + quantities::DebugString(total, 3) + u8"; ∛-2 = " +
                 quantities::DebugString(cbrt(-2)));
}

__declspec(noinline) void BenchmarkCbrtKeplerThroughput(benchmark::State& state, double (*cbrt)(double)) {
  static double zero_cycles = quantities::Infinity<double>;
  double total = 0;
  double min_cycles = quantities::Infinity<double>;
  int iterations = 0;
  constexpr std::int64_t n = 64;
  constexpr std::uint64_t low = 0x3FF0000000000000;   // 1.
  constexpr std::uint64_t high = 0x4000000000000000;  // 2.
  std::linear_congruential_engine<std::uint64_t,
                                  6364136223846793005,
                                  1442695040888963407,
                                  0>
      rng(1729);
  struct Input {
    double e;
    double E;
  };
  std::array<Input, static_cast<std::size_t>(n)> inputs;
  for (std::int64_t i = 0; i < n; ++i) {
    std::uint64_t const Y = rng() % (high - low) + low;
    double const y = principia::numerics::to_double(Y);
    std::uint64_t const Z = rng() % (high - low) + low;
    double const z = principia::numerics::to_double(Y);
    inputs[i] = {.e = 0.999 * (y - 1), .E = z - 1};
  }
  while (state.KeepRunning()) {
    auto const start = __rdtsc();
    for (std::int64_t i = 0; i < n; ++i) {
      double const e = inputs[i].e;
      double E = inputs[i].E;
      // Make sure we are in [0, 1] even if the cube root secretly the identity,
      // in constant time.
      E = to_double((to_integer(E) & ~0xFFF0000000000000) |
                    0x3FF0000000000000) - 1;
      double const s = std::sin(E);
      double const c = std::cos(E);
      // ⁰¹²³⁴⁵⁶⁷⁸⁹
      double const e² = e * e;
      double const e³ = e² * e;
      double const s² = s * s;
      double const s³ = s² * s;
      double const s⁴ = s² * s²;
      double const c² = c * c;
      double const ec = e * c;
      double const e²c² = ec * ec;
      double const R =
          cbrt(3 * e² * s * c - e³ * s³ +
               quantities::Sqrt(e³ * c² *
                                (8 * quantities::Pow<3>(1 - ec) -
                                 3 * e*s² * (1 + 4 * ec * (ec - 2)) -
                                 6 * e³ * s⁴)));
      inputs[i].E +=
          (-2 * ec + 2 * e²c² + R * R - e * s * R + e² * s²) / (ec * R);
    }
    auto const stop = __rdtsc();
    min_cycles = std::min(min_cycles, static_cast<double>(stop - start));
    ++iterations;
  }
  for (std::int64_t i = 0; i < n; ++i) {
    total += inputs[i].E / n;
  }
  if (cbrt(5) == 5) {
    zero_cycles = std::min(zero_cycles, min_cycles / n);
  }
  state.SetLabel(std::to_string(min_cycles / n - zero_cycles) +
                 " cycles " + quantities::DebugString(total, 3) + u8"; ∛-2 = " +
                 quantities::DebugString(cbrt(-2)));
}

#define CBRT_BENCHMARKS(name, f)                            \
  void BM_Latency##name(benchmark::State& state) {          \
    BenchmarkCbrtLatency(state, (f));                       \
  }                                                         \
  BENCHMARK(BM_Latency##name);                              \
  void BM_Throughput##name(benchmark::State& state) {       \
    BenchmarkCbrtThroughput(state, (f));                    \
  }                                                         \
  BENCHMARK(BM_Throughput##name);                           \
  void BM_KeplerThroughput##name(benchmark::State& state) { \
    BenchmarkCbrtKeplerThroughput(state, (f));              \
  }                                                         \
  BENCHMARK(BM_KeplerThroughput##name)

CBRT_BENCHMARKS(NoCbrt, [](double x) { return x; });
CBRT_BENCHMARKS(StdCbrt, [](double x) { return std::cbrt(x); });
CBRT_BENCHMARKS(StdSin, [](double x) { return std::sin(x)+std::cos(x); });
CBRT_BENCHMARKS(CorrectCbrt, [](double x) { return correct_cube_root(std::abs(x)).nearest_rounding; });
CBRT_BENCHMARKS(FastCorrectCbrt, &fast_correct::cbrt);
CBRT_BENCHMARKS(R3DR6Cbrt, &r3dr6::cbrt);
CBRT_BENCHMARKS(R3DR5Cbrt, &r3dr5::cbrt);
CBRT_BENCHMARKS(I3TDR5Cbrt, &i3tdr5::cbrt);
CBRT_BENCHMARKS(R5DR4FMACbrt, &r5dr4_fma::cbrt);
CBRT_BENCHMARKS(I5DR4FMACbrt, &i5dr4_fma::cbrt);

#endif

}  // namespace numerics
}  // namespace principia
