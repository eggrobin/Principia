// .\Release\x64\benchmarks.exe --benchmark_filter=MLSL
// --benchmark_repetitions=1  // NOLINT(whitespace/line_length)

#include "numerics/global_optimization.hpp"

#include <random>

#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/optimization_test_functions.hpp"

namespace principia {
namespace numerics {

using geometry::Displacement;
using geometry::Frame;
using geometry::Vector;
using quantities::Inverse;
using quantities::Length;
using quantities::si::Metre;
using testing_utilities::Branin;
using testing_utilities::GoldsteinPrice;
using testing_utilities::𝛁Branin;
using testing_utilities::𝛁GoldsteinPrice;

using World = Frame<enum class WorldTag>;

void BM_MLSLBranin(benchmark::State& state) {
  std::int64_t const points_per_round = state.range(0);
  std::int64_t const number_of_rounds = state.range(1);

  using Optimizer =
      MultiLevelSingleLinkage<double, Displacement<World>, /*dimensions=*/2>;

  auto branin = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return Branin(x₁, x₂);
  };

  auto grad_branin = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    double const g₀ = 0;
    auto const [g₁, g₂] = 𝛁Branin(x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>({0 * Metre, 2.5 * Metre, 7.5 * Metre}),
      .vertices = {
          Displacement<World>({0 * Metre, 7.5 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 7.5 * Metre}),
      }};

  auto const tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, branin, grad_branin);

  for (auto _ : state) {
    benchmark::DoNotOptimize(optimizer.FindGlobalMinima(
        points_per_round, number_of_rounds, tolerance));
  }
}

void BM_MLSLGoldsteinPrice(benchmark::State& state) {
  std::int64_t const points_per_round = state.range(0);
  std::int64_t const number_of_rounds = state.range(1);

  using Optimizer =
      MultiLevelSingleLinkage<double, Displacement<World>, /*dimensions=*/2>;

  auto goldstein_price = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    return GoldsteinPrice(x₁, x₂);
  };

  auto grad_goldstein_price = [](Displacement<World> const& displacement) {
    auto const& coordinates = displacement.coordinates();
    double const x₁ = coordinates[1] / Metre;
    double const x₂ = coordinates[2] / Metre;
    double const g₀ = 0;
    auto const [g₁, g₂] = 𝛁GoldsteinPrice(x₁, x₂);
    return Vector<Inverse<Length>, World>({g₀ / Metre, g₁ / Metre, g₂ / Metre});
  };

  Optimizer::Box const box = {
      .centre = Displacement<World>(),
      .vertices = {
          Displacement<World>({0 * Metre, 2 * Metre, 0 * Metre}),
          Displacement<World>({0 * Metre, 0 * Metre, 2 * Metre}),
      }};

  auto const tolerance = 1e-6 * Metre;
  Optimizer optimizer(box, goldstein_price, grad_goldstein_price);

  for (auto _ : state) {
    benchmark::DoNotOptimize(optimizer.FindGlobalMinima(
        points_per_round, number_of_rounds, tolerance));
  }
}

BENCHMARK(BM_MLSLBranin)->ArgsProduct({{10, 20, 50}, {10, 20, 50}});
BENCHMARK(BM_MLSLGoldsteinPrice)->ArgsProduct({{10, 20, 50}, {10, 20, 50}});

}  // namespace numerics
}  // namespace principia
