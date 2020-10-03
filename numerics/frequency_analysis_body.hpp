
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
         int wdegree_,
         template<typename, typename, int> class Evaluator>
AngularFrequency PreciseMode(
    Interval<AngularFrequency> const& fft_mode,
    Function const& function,
    PoissonSeries<double, wdegree_, Evaluator> const& weight) {
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

template<int degree_,
         typename Function,
         int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
Projection(Function const& function,
           AngularFrequency const& ω,
           PoissonSeries<double, wdegree_, Evaluator> const& weight,
           Instant const& t_min,
           Instant const& t_max) {
  std::optional<AngularFrequency> optional_ω = ω;

  // A calculator that returns optional_ω once and then stops.
  auto angular_frequency_calculator = [&optional_ω](auto const& residual) {
    auto const result = optional_ω;
    optional_ω = std::nullopt;
    return result;
  };

  return IncrementalProjection<degree_>(function,
                                        angular_frequency_calculator,
                                        weight,
                                        t_min, t_max);
}

template<int degree_,
         int periodic_degree,
         typename Function,
         typename AngularFrequencyCalculator, int wdegree_,
         template<typename, typename, int> class Evaluator>
PoissonSeries<std::invoke_result_t<Function, Instant>, degree_, Evaluator>
IncrementalProjection(Function const& function,
                      AngularFrequencyCalculator const& calculator,
                      PoissonSeries<double, wdegree_, Evaluator> const& weight,
                      Instant const& t_min,
                      Instant const& t_max,
                      mathematica::Logger* logger,
                      std::string_view celestial,
                      int interval_index) {
  using Value = std::invoke_result_t<Function, Instant>;
  using Norm = typename Hilbert<Value>::NormType;
  using Norm² = typename Hilbert<Value>::InnerProductType;
  using Normalized = typename Hilbert<Value>::NormalizedType;
  using Series = PoissonSeries<Value, degree_, Evaluator>;

  // This code follows [Kud07], section 2.  Our indices start at 0, unlike those
  // of Кудрявцев which start at 1.

  Instant const& t0 = weight.origin();

  std::optional<AngularFrequency> ω = calculator(function);
  CHECK(ω.has_value());

  std::vector<Series> basis;

  int basis_size;
  if (ω.value() == AngularFrequency{}) {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
            t0);
    basis_size = std::tuple_size_v<decltype(ω_basis)>;
    std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
  } else {
    auto const ω_basis =
        PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
            ω.value(), t0);
    constexpr int periodic_basis_size =
        2 * Hilbert<Value>::dimension * periodic_degree;
    static_assert(periodic_basis_size < std::tuple_size_v<decltype(ω_basis)>);
    basis_size = periodic_basis_size;
    std::move(ω_basis.begin(),
              ω_basis.begin() + periodic_basis_size,
              std::back_inserter(basis));
  }

  UnboundedLowerTriangularMatrix<Inverse<Norm>> α(basis_size, uninitialized);

  // Only indices 0 to m - 1 are used in this array.  At the beginning of
  // iteration m it contains Aⱼ⁽ᵐ⁻¹⁾.
  UnboundedVector<double> A(basis_size, uninitialized);

  Norm² const F₀ = InnerProduct(function, basis[0], weight, t_min, t_max);
  Norm² const Q₀₀ = InnerProduct(basis[0], basis[0], weight, t_min, t_max);
  α[0][0] = 1 / Sqrt(Q₀₀);
  A[0] = F₀ / Q₀₀;

  // At the beginning of iteration m this contains fₘ₋₁.
  auto f = function - A[0] * basis[0];

  int m_begin = 1;
  int mm = 1;
  std::vector<PoissonSeries<Normalized, degree_, Evaluator>> b;
  for (;;) {
    for (int m = m_begin; m < basis_size; ++m, ++mm) {
      // Contains Fₘ.
      Norm² const F = InnerProduct(f, basis[m], weight, t_min, t_max);

      // This vector contains Qₘⱼ.
      UnboundedVector<Norm²> Q(m + 1, uninitialized);
      for (int j = 0; j <= m; ++j) {
        Q[j] = InnerProduct(basis[m], basis[j], weight, t_min, t_max);
      }

      // This vector contains Bⱼ⁽ᵐ⁾.
      UnboundedVector<Norm> B(m, uninitialized);
      for (int j = 0; j < m; ++j) {
        Norm Σ_αⱼₛ_Qₘₛ{};
        for (int s = 0; s <= j; ++s) {
          Σ_αⱼₛ_Qₘₛ += α[j][s] * Q[s];
        }
        B[j] = -Σ_αⱼₛ_Qₘₛ;
      }

      {
        Norm² Σ_Bₛ⁽ᵐ⁾²{};
        for (int s = 0; s < m; ++s) {
          Σ_Bₛ⁽ᵐ⁾² += B[s] * B[s];
        }
        if (false && (Q[m] <= Σ_Bₛ⁽ᵐ⁾² ||
            (Q[m] - Σ_Bₛ⁽ᵐ⁾²) / std::max(Q[m], Σ_Bₛ⁽ᵐ⁾²) < 0x1.0p-24)) {
          // We arrive here when the norm of Σₛ Bₛ⁽ᵐ⁾bₛ + eₘ is small (see
          // [SN97] for the notation) and, due to rounding errors, the computed
          // value of the square of that norm ends up negative, zero, or very
          // small.  It makes no sense to have complex numbers (or infinities)
          // here because our function is real and bounded.  But even if the
          // norm could be computed but was very small, we would end up with an
          // ill-conditioned solution.  Geometrically, we are in a situation
          // where eₘ is very close to the space spanned by the (bₛ), that is,
          // by the (eₛ) for i < m.  The fact that the basis elements and no
          // longer independent when the degree increases is duely noted by
          // [CV84].  Given that eₘ effectively doesn't have benefit for the
          // projection, we just drop it and continue with the algorithm.
          LOG(ERROR) << "Q[m]: " << Q[m] << " Σ_Bₛ⁽ᵐ⁾² " << Σ_Bₛ⁽ᵐ⁾²
                     << " difference: " << Q[m] - Σ_Bₛ⁽ᵐ⁾²;
          LOG(ERROR) << "Dropping " << basis[m];
          int const basis_remaining = basis_size - m - 1;
          basis.erase(basis.begin() + m);
          α.EraseToEnd(m);
          α.Extend(basis_remaining, uninitialized);
          A.EraseToEnd(m);
          A.Extend(basis_remaining, uninitialized);
          --basis_size;
          --m;
          continue;
        } else {
          α[m][m] = 1 / Sqrt(Q[m] - Σ_Bₛ⁽ᵐ⁾²);
        }
      }

      for (int j = 0; j < m; ++j) {
        double Σ_Bₛ⁽ᵐ⁾_αₛⱼ = 0;
        for (int s = j; s < m; ++s) {
          Σ_Bₛ⁽ᵐ⁾_αₛⱼ += B[s] * α[s][j];
        }
        α[m][j] = α[m][m] * Σ_Bₛ⁽ᵐ⁾_αₛⱼ;
      }

        PoissonSeries<Normalized, degree_, Evaluator> Σ_αₘᵢ_eᵢ =
          α[m][0] * basis[0];
      for (int i = 1; i <= m; ++i) {
        Σ_αₘᵢ_eᵢ += α[m][i] * basis[i];
      }
      {
        b.push_back(Σ_αₘᵢ_eᵢ);
        for (int i = 0; i < b.size(); ++i) {
          logger->Append(absl::StrCat("bInnerProducts[",
                                      celestial,
                                      ",",
                                      interval_index,
                                      ",",
                                      mm,
                                      "]"),
                         InnerProduct(b[i], b.back(), weight, t_min, t_max));
        }
        double error = 0;
        for (int i = 0; i < b.size() - 1; ++i) {
          error += quantities::Pow<2>(
              InnerProduct(b[i], b.back(), weight, t_min, t_max));
        }
        error += quantities::Pow<2>(
            1 - InnerProduct(b.back(), b.back(), weight, t_min, t_max));
        if (Sqrt(error / m) > 0x1p-24) {
          LOG(ERROR) << "DROP " << mm;
          int const basis_remaining = basis_size - m - 1;
          basis.erase(basis.begin() + m);
          α.EraseToEnd(m);
          α.Extend(basis_remaining, uninitialized);
          A.EraseToEnd(m);
          A.Extend(basis_remaining, uninitialized);
          b.pop_back();
          --basis_size;
          --m;
          continue;
        }
      }

      A[m] = α[m][m] * α[m][m] * F;

      for (int j = 0; j < m; ++j) {
        A[j] += α[m][m] * α[m][j] * F;
      }

      {
        f -= α[m][m] * F * Σ_αₘᵢ_eᵢ;
      }
      if (logger != nullptr) {
        Norm max_residual{};
        Norm² sum_of_square_residuals{};
        for (int i = 0; i < 1 << 10; ++i) {
          auto const residual = f(t_min + i * ((t_max - t_min) / (1 << 10)));
          max_residual = std::max(max_residual, residual.Norm());
          sum_of_square_residuals += residual.Norm²();
        }
        LOG(INFO) << "max_residual=" << max_residual;
        logger->Set(
            absl::StrCat(
                "maxResidual[", celestial, ",", interval_index, ",", mm, "]"),
            max_residual,
            mathematica::ExpressIn(quantities::si::Unit<Norm>));
        logger->Set(
            absl::StrCat(
                "rmsResidual[", celestial, ",", interval_index, ",", mm, "]"),
            Sqrt(sum_of_square_residuals / (1 << 10)),
            mathematica::ExpressIn(quantities::si::Unit<Norm>));
      }
    }

    PoissonSeries<Value, degree_, Evaluator> result = A[0] * basis[0];
    for (int i = 1; i < basis_size; ++i) {
      result += A[i] * basis[i];
    }

    ω = calculator(f);
    if (!ω.has_value()) {
      return result;
    }

    int ω_basis_size;
    if (ω.value() == AngularFrequency{}) {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
              t0);
      ω_basis_size = std::tuple_size_v<decltype(ω_basis)>;
      std::move(ω_basis.begin(), ω_basis.end(), std::back_inserter(basis));
    } else {
      auto const ω_basis =
          PoissonSeriesBasisGenerator<Series, Hilbert<Value>::dimension>::Basis(
              ω.value(), t0);
      constexpr int periodic_basis_size =
          2 * Hilbert<Value>::dimension * periodic_degree;
      static_assert(periodic_basis_size < std::tuple_size_v<decltype(ω_basis)>);
      std::move(ω_basis.begin(),
                ω_basis.begin() + periodic_basis_size,
                std::back_inserter(basis));
      ω_basis_size = periodic_basis_size;
    }
    α.Extend(ω_basis_size, uninitialized);
    A.Extend(ω_basis_size, uninitialized);
    m_begin = basis_size;
    basis_size += ω_basis_size;
  }
}

}  // namespace internal_frequency_analysis
}  // namespace frequency_analysis
}  // namespace numerics
}  // namespace principia
