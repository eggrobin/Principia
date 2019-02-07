#pragma once

#include "quantities/uncertainty.hpp"

namespace principia {
namespace quantities {
namespace internal_uncertainty {

std::ostream& operator<<(std::ostream& out,
                         MeasurementResult<double> measurement) {
  // REMOVE BEFORE FLIGHT: This needs to not summon UB on infinities or NaNs.
  int const value_decimal_exponent =
      std::floor(std::log10(std::abs(measurement.measured_value)));
  int const uncertainty_decimal_exponent =
      std::floor(std::log10(measurement.standard_uncertainty));
  int const significant_digits =
      value_decimal_exponent - uncertainty_decimal_exponent;
  std::string value_digits = std::to_string(static_cast<int>(std::nearbyint(
      std::abs(measurement.measured_value) *
      std::pow(10, (significant_digits + 2) - value_decimal_exponent - 1))));
  std::string uncertainty_digits = std::to_string(static_cast<int>(
      std::nearbyint(std::abs(measurement.standard_uncertainty) *
                     std::pow(10, 2 - uncertainty_decimal_exponent - 1))));
  if (value_decimal_exponent >= 0 &&
      significant_digits >= value_decimal_exponent + 1) {
    return out << (std::signbit(measurement.measured_value) ? "-" : "+")
               << value_digits.substr(0, value_decimal_exponent + 1) << "."
               << value_digits.substr(value_decimal_exponent + 1) << "("
               << uncertainty_digits << ")";
  }
  return out << (std::signbit(measurement.measured_value) ? "-" : "+")
             << value_digits[0] << "." << value_digits.substr(1) << "("
             << uncertainty_digits << u8") × 10^"
             << (value_decimal_exponent >= 0 ? "+" : "")
             << value_decimal_exponent;
}

template<typename T, typename U>
MeasurementResult<Sum<T, U>> operator+(MeasurementResult<T> const& left,
                                       U const& right) {
  return {left.measured_value + right, left.standard_uncertainty};
}
template<typename T, typename U>
MeasurementResult<Sum<T, U>> operator+(T const& left,
                                       MeasurementResult<U> const& right) {
  return {left + right.measured_value, right.standard_uncertainty};
}

template<typename T, typename U>
MeasurementResult<Difference<T, U>> operator-(MeasurementResult<T> const& left,
                                              U const& right) {
  return {left.measured_value - right, left.standard_uncertainty};
}

template<typename T, typename U>
MeasurementResult<Difference<T, U>> operator-(
    T const& left,
    MeasurementResult<U> const& right) {
  return {left - right.measured_value, right.standard_uncertainty};
}

template<typename T, typename U>
MeasurementResult<Product<T, U>> operator*(MeasurementResult<T> const& left,
                                           U const& right) {
  return {left.measured_value * right,
          left.standard_uncertainty * Abs(right)};
}

template<typename T, typename U>
MeasurementResult<Product<T, U>> operator*(T const& left,
                                           MeasurementResult<U> const& right) {
  return {left * right.measured_value,
          Abs(left) * right.standard_uncertainty};
}

template<typename T, typename U>
MeasurementResult<Quotient<T, U>> operator/(MeasurementResult<T> const& left,
                                            U const& right) {
  return {left.measured_value / right,
          left.standard_uncertainty / Abs(right)};
}

template<typename T, typename U>
MeasurementResult<Quotient<T, U>> operator/(T const& left,
                                            MeasurementResult<U> const& right) {
  return {
      left / right.measured_value,
      Abs(left) * right.standard_uncertainty / Pow<2>(right.measured_value)};
}

template<typename T, typename U>
MeasurementResult<Quotient<T, U>> operator/(MeasurementResult<T> const& left,
                                            MeasurementResult<U> const& right) {
  return {
      left.measured_value / right.measured_value,
      Sqrt(Pow<2>(left.standard_uncertainty) / Pow<2>(right.measured_value) +
           Pow<2>(right.standard_uncertainty) * Pow<2>(left.measured_value) /
               Pow<4>(left.measured_value))};
}

template<typename T>
T Average(std::vector<T> const& samples) {
  // TODO(egg): we should have a version of Barycentre without coefficients.
  Difference<T> total{};
  for (T const& sample : samples) {
    total += sample - T{};
  }
  return total / samples.size() + T{};
}

template<typename T>
MeasurementResult<T> AverageOfIndependent(std::vector<T> const& measured_values) {
  // We use the notation from GUM 4.2.
  int const n = time_series.size();
  auto const& q = measured_values;
  T const q̄ = Average(q);
  // The estimate of the variance.
  Square<Difference<T>> s²_q{};
  for (T const& qⱼ : q) {
    s²_q += Pow<2>(qⱼ - q̄);
  }
  s²_q /= n - 1;
  // The experimental standard deviation of the mean.
  Difference<T> const s_q̄ = Sqrt(s²_q / n);
  return {q̄, s_q̄};
}

template<typename T>
MeasurementResult<T> AverageOfCorrelated(std::vector<T> const& time_series) {
  // See:
  // — Grant Foster (1996), Time Series Analysis by Projection. I. Statistical
  //   Properties of Fourier Analysis, equation 5.4;
  // — Chris Koen and Fred Lombard (1993), The analysis of indexed astronomical
  //   time series — I. Basic methods, equations 15, 19, and 2.

  int const n = time_series.size();
  auto const& y = time_series;
  T const ȳ = Average(y);

  // TODO(egg): L = 1 works well when the time series is dominated by periodic
  // signals rather than noise; for a noisier time series, higher values of L
  // behave better; perhaps to be safe we should use a slightly higher L.
  // For L = n - 1, the estimate of the spectral density at 0 frequency becomes
  // the experimental variance.
  constexpr int L = 1;

  // The estimated spectral density at zero frequency, see Koen and Lombard
  // (1993), section 4.4.
  Square<Difference<T>> Ŝε_0{};
  for (int j = 1; j <= L; ++j) {
    // The real and imaginary parts of √n F(y)ⱼ, where y is the time series and F
    // is the unitary discrete Fourier transform.
    Difference<T> sqrt_n_re_fourier_j{};
    Difference<T> sqrt_n_im_fourier_j{};
    // We use the convention from Koen and Lombard (2004), ω = 2πj/n, rather
    // than ω = j/n as in Koen and Lombard (1993).
    double const ω = 2 * π * j / n;
    for (int t = 0; t < n; ++t) {
      // Note that there is a typo in equation 2 of Koen & Lombard, the t is
      // missing from the exponent.  A correct expression for the periodogram
      // may be found in equation 4 of Chris Koen and Fred Lombard (2004), The
      // analysis of indexed astronomical time series — IX. A period change
      // test.
      Difference<T> const Δy = y[t] - ȳ;
      sqrt_n_re_fourier_j += Δy * std::cos(ω * t);
      sqrt_n_im_fourier_j += Δy * std::sin(ω * t);
    }
    Ŝε_0 += Pow<2>(sqrt_n_re_fourier_j) + Pow<2>(sqrt_n_im_fourier_j);
  }
  Ŝε_0 /= (n * L);

  return {ȳ, Sqrt(Ŝε_0 / n)};
}

template<typename Argument, typename Value>
LinearModel<Argument, Value> LinearRegression(
    std::vector<Argument> const& arguments,
    std::vector<Value> const& values) {
  CHECK_EQ(arguments.size(), values.size());
  int const n = arguments.size();
  // We follow the notation from Foster (1996).
  auto const& t = arguments;
  auto const& x = values;
  // ⟨𝟏|t⟩.
  auto const t̄ = Average(t);
  // ⟨𝟏|x⟩ = ⟨𝟏|x⟩ / ⟨𝟏|𝟏⟩, the projection onto the constant trial function.
  auto const x̄ = Average(x);

  // The variation Vt of t is the inner product of
  // |t⟩ - ⟨𝟏|t⟩|𝟏⟩ = |t⟩ - t̄|𝟏⟩ ≕ |l⟩ with itself.
  // If |Argument| is a vector space rather than an affine space,
  // this simplifies to ⟨t|t⟩ - ⟨𝟏|t⟩², see equation 2.11.  Note that this is
  // not quite the estimated variance of t, which is n Vt / (n-1), see
  // equation 2.12.
  Square<Difference<Argument>> variation_t{};
  for (int α = 0; α < n; ++α) {
    variation_t += Pow<2>(t[α] - t̄);
  }
  variation_t /= n;
  // ⟨l|x⟩.
  Product<Difference<Value>, Difference<Argument>> inner_product{};
  for (int α = 0; α < n; ++α) {
    inner_product += (t[α] - t̄) *  (x[α] - x̄);
  }
  inner_product /= n;
  // ⟨l|x⟩ / ⟨l|l⟩, the projection of |x⟩ onto the linear trial function
  // |l⟩, i.e., the slope.
  Derivative<Value, Argument> const y = inner_product / variation_t;

  LinearModel<Argument, Value> result;
  result.mean_value.measured_value = x̄;
  result.slope.measured_value = y;

  // Now we evaluate the uncertainties.

  // s² under the strong signal hypothesis, from equation 6.4.  This corresponds
  // to GUM H.13f.
  Square<Difference<Value>> s²{};
  for (int α = 0; α < n; ++α) {
    Difference<Value> const residual = x[α] - (y * (t[α] - t̄) + x̄);
    s² += Pow<2>(residual);
  }
  s² /= n - 2;

  result.mean_value.standard_uncertainty = Sqrt(s² / n);
  result.slope.standard_uncertainty = Sqrt(s² / (n * variation_t));
  return result;
}

}  // namespace internal_uncertainty
}  // namespace quantities
}  // namespace principia
