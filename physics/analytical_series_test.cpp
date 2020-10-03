
#include <algorithm>
#include <chrono>
#include <vector>

#include "absl/strings/str_cat.h"
#include "astronomy/frames.hpp"
#include "geometry/interval.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/apodization.hpp"
#include "numerics/fast_fourier_transform.hpp"
#include "numerics/frequency_analysis.hpp"
#include "numerics/poisson_series.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"

namespace principia {
namespace physics {

using astronomy::ICRS;
using geometry::Displacement;
using geometry::Instant;
using geometry::Interval;
using geometry::Position;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::QuinlanTremaine1990Order12;
using numerics::EstrinEvaluator;
using numerics::FastFourierTransform;
using numerics::PoissonSeries;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Length;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;

namespace apodization = numerics::apodization;
namespace frequency_analysis = numerics::frequency_analysis;

static constexpr int secular_degree = 5;
static constexpr int periodic_degree = 2;
static constexpr int log2_number_of_samples = 10;
static constexpr int number_of_frequencies = 20;

class AnalyticalSeriesTest : public ::testing::Test {
 protected:
  AnalyticalSeriesTest()
      : logger_(TEMP_DIR / "analytical_series.wl",
               /*make_unique=*/false) {
    google::LogToStderr();
    logger_.Set("secularDegree", secular_degree);
    logger_.Set("periodicDegree", periodic_degree);
    logger_.Set("frequencies", number_of_frequencies);
  }

  template<int degree>
  std::unique_ptr<
      PoissonSeries<Displacement<ICRS>, secular_degree, EstrinEvaluator>>
  ComputeCompactRepresentation(std::string_view const celestial,
                               ContinuousTrajectory<ICRS> const& trajectory) {
    Instant const t_min = trajectory.t_min();
    std::unique_ptr<PoissonSeries<Displacement<ICRS>,
                                  secular_degree,
                                  EstrinEvaluator>>
        result;
    int interval_index = 1;
    for (Instant t_max = t_min + 80 * Minute; t_max < trajectory.t_max();
         t_max = t_min + (t_max - t_min) * 2, ++interval_index) {
    auto const piecewise_poisson_series =
        trajectory.ToPiecewisePoissonSeries<degree>(t_min, t_max);

    int step = 0;

    auto angular_frequency_calculator =
        [this, celestial, &step, t_min, t_max](
            auto const& residual) -> std::optional<AngularFrequency> {
      Time const Δt = (t_max - t_min) / (1 << log2_number_of_samples);
      if (step == 0) {
        ++step;
        return AngularFrequency();
      } else if (step <= number_of_frequencies) {
        ++step;
        Length max_residual;
        std::vector<Displacement<ICRS>> residuals;
        for (int i = 0; i < 1 << log2_number_of_samples; ++i) {
          residuals.push_back(residual(t_min + i * Δt));
          max_residual = std::max(max_residual, residuals.back().Norm());
        }
        auto fft =
            std::make_unique<FastFourierTransform<Displacement<ICRS>,
                                                  1 << log2_number_of_samples>>(
                residuals, Δt);
        auto const mode = fft->Mode();
        Interval<Time> const period{2 * π * Radian / mode.max,
                                    2 * π * Radian / mode.min};
        auto const precise_mode = frequency_analysis::PreciseMode(
            mode, residual, apodization::Hann<EstrinEvaluator>(t_min, t_max));
        auto const precise_period = 2 * π * Radian / precise_mode;
        LOG(INFO) << "precise_period=" << precise_period;
        return precise_mode;
      } else {
        Length max_residual;
        for (int i = 0; i < 1 << log2_number_of_samples; ++i) {
          max_residual =
              std::max(max_residual, residual(t_min + i * Δt).Norm());
        }
        return std::nullopt;
      }
    };

    result = std::make_unique<PoissonSeries<Displacement<ICRS>,
                                            secular_degree,
                                            EstrinEvaluator>>(
        frequency_analysis::IncrementalProjection<secular_degree,
                                                  periodic_degree>(
            piecewise_poisson_series,
            angular_frequency_calculator,
            apodization::Dirichlet<EstrinEvaluator>(t_min, t_max),
            t_min,
            t_max,
            &logger_,
            celestial,
            interval_index)); 
  }
  logger_.Set("intervals", interval_index - 1);
  return result;
  }

  mathematica::Logger logger_;
};

#define PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(                   \
    degree, approximation, celestial, trajectory)                        \
  case degree: {                                                         \
    approximation = ComputeCompactRepresentation<(degree)>(celestial, trajectory);  \
    break;                                                               \
  }

#if !_DEBUG
TEST_F(AnalyticalSeriesTest, CompactRepresentation) {
  SolarSystem<ICRS> solar_system_at_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");

  // NOTE(phl): Keep these parameters aligned with
  // sol_numerics_blueprint.proto.txt.
  auto const ephemeris = solar_system_at_j2000.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                             Position<ICRS>>(),
          /*step=*/10 * Minute));
  ephemeris->Prolong(solar_system_at_j2000.epoch() + 1 * JulianYear);

  for (auto const& celestial : {"Io", "Jupiter", "Phobos", "Pluto", "Earth", "Moon", "Venus"}) {
    auto const start = std::chrono::system_clock::now();
    auto const& trajectory =
        solar_system_at_j2000.trajectory(*ephemeris, celestial);
    int const piecewise_poisson_series_degree =
        trajectory.PiecewisePoissonSeriesDegree(trajectory.t_min(),
                                                trajectory.t_max());
    std::unique_ptr<PoissonSeries<Displacement<ICRS>,
                                  secular_degree,
                                  EstrinEvaluator>> approximation;

    LOG(INFO) << "---- " << celestial << " degree "
              << piecewise_poisson_series_degree;

    switch (piecewise_poisson_series_degree) {
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          3, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          4, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          5, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          6, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          7, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          8, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          9, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          10, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          11, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          12, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          13, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          14, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          15, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          16, approximation, celestial, trajectory);
      PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE(
          17, approximation, celestial, trajectory);
      default:
        LOG(FATAL) << "Unexpected degree "
                   << piecewise_poisson_series_degree;
    };
    auto const end = std::chrono::system_clock::now();
    LOG(INFO) << "Timing " << celestial << ": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       start)
                     .count()
              << " ms";

    logger_.Set(absl::StrCat("approximation", celestial),
                *approximation,
                mathematica::ExpressIn(Metre, Second, Radian));
  }
}
#endif

#undef PRINCIPIA_COMPUTE_COMPACT_REPRESENTATION_CASE

}  // namespace physics
}  // namespace principia
