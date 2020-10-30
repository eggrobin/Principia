﻿
#pragma once

#include "numerics/frequency_analysis.hpp"

#include <algorithm>
#include <functional>
#include <vector>

#include "base/tags.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "numerics/poisson_series_basis.hpp"
#include "numerics/root_finders.hpp"
#include "numerics/unbounded_arrays.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace frequency_analysis {
namespace internal_frequency_analysis {

using base::uninitialized;
using geometry::Hilbert;
using geometry::Vector;
using quantities::Inverse;
using quantities::Sqrt;
using quantities::Square;
using quantities::SquareRoot;

template<typename Function,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double,
                  aperiodic_wdegree, periodic_wdegree,
                  Evaluator> const& weight) {
  auto const weighted_function = weight * function;
  auto const weighted_function_spectrum = weighted_function.FourierTransform();

  auto power =
      [&weighted_function_spectrum](AngularFrequency const& ω) {
        return weighted_function_spectrum(ω).Norm²();
      };

  return Brent(power,
               fft_mode.min,
               fft_mode.max,
               std::greater<>());
}

template<int aperiodic_degree, int periodic_degree,
         typename Function,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>,
              aperiodic_degree, periodic_degree,
              Evaluator>
Projection(Function const& function,
           AngularFrequency const& ω,
           PoissonSeries<double,
                         aperiodic_wdegree, periodic_wdegree,
                         Evaluator> const& weight,
           Instant const& t_min,
           Instant const& t_max) {
  std::optional<AngularFrequency> optional_ω = ω;

  // A calculator that returns optional_ω once and then stops.
  auto angular_frequency_calculator = [&optional_ω](auto const& residual) {
    auto const result = optional_ω;
    optional_ω = std::nullopt;
    return result;
  };

  return IncrementalProjection<aperiodic_degree, periodic_degree>(
      function,
      angular_frequency_calculator,
      weight,
      t_min, t_max);
}

template<int aperiodic_degree, int periodic_degree,
         typename Function,
         typename AngularFrequencyCalculator,
         int aperiodic_wdegree, int periodic_wdegree,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>,
              aperiodic_degree, periodic_degree,
              Evaluator>
IncrementalProjection(Function const& function,
                      AngularFrequencyCalculator const& calculator,
                      PoissonSeries<double,
                                    aperiodic_wdegree, periodic_wdegree,
                                    Evaluator> const& weight,
                      Instant const& t_min,
                      Instant const& t_max,
                      std::string const& tag,
                      mathematica::Logger& logger) {
  using Value = std::invoke_result_t<Function, Instant>;
  using Norm = typename Hilbert<Value>::NormType;
  using Norm² = typename Hilbert<Value>::Norm²Type;
  using Normalized = typename Hilbert<Value>::NormalizedType;
  using Series = PoissonSeries<Value,
                               aperiodic_degree, periodic_degree,
                               Evaluator>;

  // This code follows [Kud07], section 2.  Our indices start at 0, unlike those
  // of Кудрявцев which start at 1.

  Instant const& t0 = weight.origin();

  std::optional<AngularFrequency> ω = calculator(function);
  CHECK(ω.has_value());

  std::vector<Series> basis;

  int basis_size;
  if (ω.value() == AngularFrequency{}) {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series,
                                    Hilbert<Value>::dimension,
                                    aperiodic_degree>::Basis(t0);
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
  } else {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series,
                                    Hilbert<Value>::dimension,
                                    periodic_degree>::Basis(ω.value(), t0);
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
  }

  // This is logically Q in the QR decomposition of basis.
  std::vector<PoissonSeries<Normalized,
                            aperiodic_degree, periodic_degree,
                            Evaluator>> q;

  auto const a₀ = basis[0];
  auto const r₀₀ = a₀.Norm(weight, t_min, t_max);
  q.push_back(a₀ / r₀₀);

  auto const A₀ = InnerProduct(function, q[0], weight, t_min, t_max);

  Series F = A₀ * q[0];
  auto f = function - F;

  int m_begin = 1;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m) {
      auto aₘ⁽ᵏ⁾ = basis[m];
      for (int k = 0; k < m; ++k) {
        auto const rₖₘ = InnerProduct(q[k], aₘ⁽ᵏ⁾, weight, t_min, t_max);
        aₘ⁽ᵏ⁾ -= rₖₘ * q[k];
      }

      auto const rₘₘ = aₘ⁽ᵏ⁾.Norm(weight, t_min, t_max);
      q.push_back(aₘ⁽ᵏ⁾ / rₘₘ);
      DCHECK_EQ(m + 1, q.size());

      Norm const Aₘ = InnerProduct(f, q[m], weight, t_min, t_max);

      f -= Aₘ * q[m];
      F += Aₘ * q[m];
      Norm max_residual;
      Norm² sum_of_squared_residuals;
static constexpr int log2_number_of_samples = 10;
auto const Δt = (t_max - t_min) / (1 << log2_number_of_samples);
      for (int i = 0; i < 1 << log2_number_of_samples; ++i) {
        max_residual =
            std::max(max_residual, Hilbert<Value>::Norm(f(t_min + i * Δt)));
        sum_of_squared_residuals += Hilbert<Value>::Norm²(f(t_min + i * Δt));
      }
      logger.Set(absl::StrCat("maxResidual[", tag, ",", m, "]"),
                 max_residual,
                 mathematica::ExpressIn(quantities::si::Metre));
      logger.Set(absl::StrCat("rmsResidual[", tag, ",", m, "]"),
                 quantities::Sqrt(sum_of_squared_residuals /
                                  (1 << log2_number_of_samples)),
                 mathematica::ExpressIn(quantities::si::Metre));
    }

    ω = calculator(f);
    if (!ω.has_value()) {
      return F;
    }

    int ω_basis_size;
    if (ω.value() == AngularFrequency{}) {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series,
                                      Hilbert<Value>::dimension,
                                      aperiodic_degree>::Basis(t0);
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    } else {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series,
                                      Hilbert<Value>::dimension,
                                      periodic_degree>::Basis(ω.value(), t0);
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    }
    m_begin = basis_size;
    basis_size += ω_basis_size;
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
