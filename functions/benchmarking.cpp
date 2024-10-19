#include "functions/benchmarking.hpp"

#include "absl/strings/str_format.h"
#include "geometry/sign.hpp"
#include "numerics/root_finders.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace functions {
namespace _benchmarking {

using namespace principia::geometry::_sign;
using namespace principia::numerics::_root_finders;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_quantities;

  std::string MeasurementResult::ToGUMString() const {
  if (standard_uncertainty == 0) {
    return DebugString(value);
  }
  double const floor_log10_u = std::floor(std::log10(standard_uncertainty));
  std::int64_t value_integer_digits = std::floor(std::log10(value)) + 1;
  std::int64_t uncertainty_digits = 2;
  std::int64_t digits_shown =
      std::floor(std::log10(value)) - floor_log10_u + uncertainty_digits;
  std::int64_t fractional_digits_shown = digits_shown - value_integer_digits;
  if (fractional_digits_shown < 0) {
    digits_shown += -fractional_digits_shown;
    uncertainty_digits += -fractional_digits_shown;
    fractional_digits_shown = 0;
  }
  if (fractional_digits_shown == 1) {
    CHECK_EQ(uncertainty_digits, 2);
    double uncertainty_parenthetical =
        std::ceil(10 * standard_uncertainty) / 10;
    return absl::StrFormat("%.1f(%03.1f)", value, uncertainty_parenthetical);
  } else {
    std::int64_t uncertainty_parenthetical =
        std::ceil(standard_uncertainty *
                  std::pow(10, uncertainty_digits - 1 - floor_log10_u));
    return absl::StrFormat("%.*f(%0*d)",
                           fractional_digits_shown,
                           value,
                           uncertainty_digits,
                           uncertainty_parenthetical);
  }
}

// From [Coh51].
MeasurementResult LogNormalTerminus(std::vector<double> const& x) {
  if (x.empty()) {
    return {0, 0};
  }
  if (x.size() == 1) {
    return {x[0], 0};
  }
  double const n = x.size();
  double const n² = n * n;
  auto const Σ₁ⁿ = [n](auto const summand) {
    double Σ = 0;
    for (int i = 0; i < n; ++i) {
      Σ += summand(i);
    }
    return Σ;
  };
  auto const λ = [&Σ₁ⁿ, &x, n, n²](double α) {
    double const Σ₁ⁿlog_xᵢ_minus_α =
        Σ₁ⁿ([&](int i) { return std::log(x[i] - α); });
    double const Σ₁ⁿlog²_xᵢ_minus_α =
        Σ₁ⁿ([&](int i) { return Pow<2>(std::log(x[i] - α)); });
    return Σ₁ⁿ([&](int i) { return 1 / (x[i] - α); }) *
               (n * Σ₁ⁿlog_xᵢ_minus_α - n * Σ₁ⁿlog²_xᵢ_minus_α +
                Pow<2>(Σ₁ⁿlog_xᵢ_minus_α)) -
           n² * Σ₁ⁿ([&](int i) { return std::log(x[i] - α) / (x[i] - α); });
  };
  Sign sign_λ_0(λ(0));
  // λ(x₁) is NaN, and λ is so ill-conditioned there that it has the wrong
  // sign just below, so we use a cheesy factor.
  double const x₁ = *std::min_element(x.begin(), x.end());
  double cheese = 1;
  while (Sign(λ((1 - cheese) * x₁)) == sign_λ_0) {
    cheese /= 2;
    if (cheese < 0x1p-53) {
      // The MLE is very close to the minimum; in that limit the variance
      // becomes 0.
      return {.value = x₁, .standard_uncertainty = 0};
    }
  }
  double const α = Brent(λ, 0.0, x₁ * (1 - cheese));
  double const Σ₁ⁿlog_xᵢ_minus_α =
      Σ₁ⁿ([&](int i) { return std::log(x[i] - α); });
  double const Σ₁ⁿlog²_xᵢ_minus_α =
      Σ₁ⁿ([&](int i) { return Pow<2>(std::log(x[i] - α)); });
  double const β = std::exp(1 / n * Σ₁ⁿlog_xᵢ_minus_α);
  double const γ² =
      1 / n * Σ₁ⁿlog²_xᵢ_minus_α - Pow<2>(1 / n * Σ₁ⁿlog_xᵢ_minus_α);
  double const ω = std::exp(γ²);
  double const β² = β * β;
  double const α_variance = β² * γ² / (n * ω * (ω * (1 + γ²) - 2 * γ² - 1));
  return {.value = α, .standard_uncertainty = Sqrt(α_variance)};
}

__declspec(noinline) double __cdecl identity(double x) {
  return x;
}

__declspec(noinline) MeasurementResult
    BenchmarkFunctionThroughput(
    double (__cdecl *f)(double),
    std::function<double()> get_input,
    std::int64_t const samples,
    MeasurementResult const identity_throughput) {
  std::vector<double> cycle_counts;
  constexpr std::int64_t n = 1 << 16;
  std::array<double, static_cast<std::size_t>(n)> inputs{};
  for (std::int64_t j = 0; j < samples; ++j) {
    for (std::int64_t i = 0; i < n; ++i) {
      inputs[i] = (inputs[i] + get_input()) - inputs[i];
    }
    auto const start = __rdtsc();
    for (std::int64_t i = 0; i < n; ++i) {
      double const result = f(inputs[i]);
      inputs[i] = f(inputs[i]);
    }
    auto const stop = __rdtsc();
    cycle_counts.push_back((double)(stop - start) / n);
  }
  MeasurementResult throughput = LogNormalTerminus(cycle_counts);
  throughput = {.value = throughput.value - identity_throughput.value,
                .standard_uncertainty =
                    Sqrt(Pow<2>(identity_throughput.standard_uncertainty) +
                         Pow<2>(identity_throughput.standard_uncertainty))};
  return throughput;
}

MeasurementResult BenchmarkFunctionThroughput(double (__cdecl *f)(double),
                                              std::function<double()> get_input,
                                              std::int64_t const samples) {
  return BenchmarkFunctionThroughput(
      f,
      get_input,
      samples,
      BenchmarkFunctionThroughput(
          &identity, get_input, samples, MeasurementResult{0, 0}));
}

__declspec(noinline) MeasurementResult
    BenchmarkFunctionLatency(
    double (__cdecl *f)(double),
    std::function<double()> get_input,
    std::int64_t const samples,
    MeasurementResult const identity_latency) {
  std::vector<double> cycle_counts;
  constexpr std::int64_t n = 1 << 16;
  for (int j = 0; j < samples; ++j) {
    std::array<double, static_cast<std::size_t>(n)> inputs;
    for (std::int64_t i = 0; i < n; ++i) {
      inputs[i] = get_input();
    }
    auto const start = __rdtsc();
    double x = inputs[0];
    for (std::int64_t i = 0; i < n; ++i) {
      double const result = f(x);
      x = result + inputs[i] - result;
    }
    auto const stop = __rdtsc();
    LOG(INFO) << x;
    cycle_counts.push_back((double)(stop - start) / n);
  }
  MeasurementResult latency = LogNormalTerminus(cycle_counts);
  latency = {.value = latency.value - identity_latency.value,
             .standard_uncertainty =
                 Sqrt(Pow<2>(latency.standard_uncertainty) +
                      Pow<2>(identity_latency.standard_uncertainty))};
  if (f == &identity) {
    LOG(ERROR) << "Identity latency:" << latency.ToGUMString();
  }
  return latency;
}

MeasurementResult BenchmarkFunctionLatency(double (__cdecl *f)(double),
                                           std::function<double()> get_input,
                                           std::int64_t const samples) {
  LOG(ERROR) << "Latency including overhead:"
             << BenchmarkFunctionLatency(
                    f, get_input, samples, MeasurementResult{0, 0}).ToGUMString();
  return BenchmarkFunctionLatency(
      f,
      get_input,
      samples,
      BenchmarkFunctionLatency(
          &identity, get_input, samples, MeasurementResult{0, 0}));
}

}  // namespace _benchmarking
}  // namespace functions
}  // namespace principia
