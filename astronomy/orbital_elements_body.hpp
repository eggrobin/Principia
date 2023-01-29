#pragma once

#include "astronomy/orbital_elements.hpp"

#include <algorithm>
#include <vector>

#include "absl/strings/str_cat.h"
#include "base/jthread.hpp"
#include "base/status_utilities.hpp"
#include "numerics/quadrature.hpp"
#include "integrators/explicit_runge_kutta_integrator.hpp"
#include "integrators/methods.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "quantities/elementary_functions.hpp"

#include "mathematica/logger.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbital_elements {

using base::this_stoppable_thread;
using geometry::Velocity;
using integrators::AdaptiveStepSizeIntegrator;
using integrators::EmbeddedExplicitRungeKuttaIntegrator;
using integrators::ExplicitLinearMultistepIntegrator;
using integrators::ExplicitRungeKuttaIntegrator;
using integrators::ExplicitFirstOrderOrdinaryDifferentialEquation;
using integrators::InitialValueProblem;
using integrators::methods::AdamsBashforthOrder4;
using integrators::methods::AdamsBashforthOrder5;
using integrators::methods::AdamsBashforthOrder6;
using integrators::methods::DormandPrince1986RK547FC;
using integrators::methods::Kutta1901Vσ1;
using numerics::quadrature::AutomaticClenshawCurtis;
using physics::DegreesOfFreedom;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using quantities::ArcTan;
using quantities::Cos;
using quantities::Mod;
using quantities::Pow;
using quantities::Product;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Square;
using quantities::Tan;
using quantities::UnwindFrom;
using quantities::si::Radian;

template<typename PrimaryCentred>
absl::StatusOr<OrbitalElements> OrbitalElements::ForTrajectory(
    Trajectory<PrimaryCentred> const& trajectory,
    MassiveBody const& primary,
    Body const& secondary,
    bool const sample_dense_elements) {
  OrbitalElements orbital_elements;
  LOG(ERROR)
      << "DiscreteTrajectory has "
      << dynamic_cast<physics::DiscreteTrajectory<PrimaryCentred> const&>(
             trajectory)
             .size()
      << " points";
  if (trajectory.t_min() >= trajectory.t_max()) {
    return absl::InvalidArgumentError(
        absl::StrCat("trajectory has min time ",
                     DebugString(trajectory.t_min()),
                     " and max time ",
                     DebugString(trajectory.t_max())));
  }

  auto osculating_elements =
      [&primary, &secondary, &trajectory](
          Instant const& time) -> KeplerianElements<PrimaryCentred> {
    DegreesOfFreedom<PrimaryCentred> const primary_dof{
        PrimaryCentred::origin, PrimaryCentred::unmoving};
    DegreesOfFreedom<PrimaryCentred> const degrees_of_freedom =
        trajectory.EvaluateDegreesOfFreedom(time);
    return KeplerOrbit<PrimaryCentred>(
               primary, secondary, degrees_of_freedom - primary_dof, time)
        .elements_at_epoch();
  };

  auto osculating_wound_λ =
      [&osculating_elements](Instant const& time) -> Angle {
    auto const elements = osculating_elements(time);
    Angle const& ϖ = *elements.longitude_of_periapsis;
    Angle const& M = *elements.mean_anomaly;
    return ϖ + M;
  };

  orbital_elements.radial_distances_ = RadialDistances(trajectory);


  KeplerianElements<PrimaryCentred> const initial_osculating_elements =
      osculating_elements(trajectory.t_min());
  Time const estimated_period = *initial_osculating_elements.period;

  std::vector<Angle> unwound_λs;
  for (Instant t = trajectory.t_min(); t <= trajectory.t_max();
       t += estimated_period / 3) {
    Angle const λ = osculating_wound_λ(t);
    unwound_λs.push_back(unwound_λs.empty() ? λ
                                            : UnwindFrom(unwound_λs.back(), λ));
  }

  auto osculating_equinoctial_elements =
      [&estimated_period,
       &osculating_elements,
       t_min = trajectory.t_min(),
       &unwound_λs](Instant const& time) -> EquinoctialElements {
    auto const elements = osculating_elements(time);
    double const& e = *elements.eccentricity;
    Angle const& ϖ = *elements.longitude_of_periapsis;
    Angle const& Ω = elements.longitude_of_ascending_node;
    Angle const& M = *elements.mean_anomaly;
    Angle const& i = elements.inclination;
    double const tg_½i = Tan(i / 2);
    double const cotg_½i = 1 / tg_½i;
    return {.t = time,
            .a = *elements.semimajor_axis,
            .h = e * Sin(ϖ),
            .k = e * Cos(ϖ),
            .λ = UnwindFrom(unwound_λs[(time - t_min) / (estimated_period / 3)],
                            ϖ + M),
            .p = tg_½i * Sin(Ω),
            .q = tg_½i * Cos(Ω),
            .pʹ = cotg_½i * Sin(Ω),
            .qʹ = cotg_½i * Cos(Ω)};
  };

  auto const sidereal_period = SiderealPeriod(
      osculating_equinoctial_elements, trajectory.t_min(), trajectory.t_max());
  RETURN_IF_ERROR(sidereal_period);
  orbital_elements.sidereal_period_ = sidereal_period.value();

  if (!IsFinite(orbital_elements.sidereal_period_) ||
      orbital_elements.sidereal_period_ <= Time{}) {
    // Guard against NaN sidereal periods (from hyperbolic orbits) or negative
    // sidereal periods (from aberrant trajectories, see #2811).
    return absl::OutOfRangeError(
        "sidereal period is " + DebugString(orbital_elements.sidereal_period_));
  }
  auto mean_equinoctial_elements =
      MeanEquinoctialElements(osculating_equinoctial_elements,
                              trajectory.t_min(),
                              trajectory.t_max(),
                              orbital_elements.sidereal_period_);
  RETURN_IF_ERROR(mean_equinoctial_elements);
  orbital_elements.mean_equinoctial_elements_ =
      std::move(mean_equinoctial_elements).value();


  if (sample_dense_elements) {
    orbital_elements.osculating_equinoctial_elements_;
    for (Instant t = trajectory.t_min(); t <= trajectory.t_max();
         t += orbital_elements.sidereal_period_ / 256) {
      orbital_elements.osculating_equinoctial_elements_.push_back(
          osculating_equinoctial_elements(t));
    }
  }

  if (orbital_elements.mean_equinoctial_elements_.size() < 2) {
    return absl::OutOfRangeError(
        "trajectory does not span one sidereal period: sidereal period is " +
            DebugString(orbital_elements.sidereal_period_) +
            ", trajectory spans " +
            DebugString(trajectory.t_max() - trajectory.t_min()));
  }
  auto mean_classical_elements =
      ToClassicalElements(orbital_elements.mean_equinoctial_elements_);
  RETURN_IF_ERROR(mean_classical_elements);
  orbital_elements.mean_classical_elements_ =
      std::move(mean_classical_elements).value();
  RETURN_IF_ERROR(orbital_elements.ComputePeriodsAndPrecession());
  RETURN_IF_ERROR(orbital_elements.ComputeIntervals());
  return orbital_elements;
}

inline std::vector<OrbitalElements::ClassicalElements> const&
OrbitalElements::mean_elements() const {
  return mean_classical_elements_;
}

inline Time OrbitalElements::sidereal_period() const {
  return sidereal_period_;
}

inline Time OrbitalElements::nodal_period() const {
  return nodal_period_;
}

inline Time OrbitalElements::anomalistic_period() const {
  return anomalistic_period_;
}

inline AngularFrequency OrbitalElements::nodal_precession() const {
  return nodal_precession_;
}

inline Interval<Length> OrbitalElements::mean_semimajor_axis_interval() const {
  return mean_semimajor_axis_interval_;
}

inline Interval<double> OrbitalElements::mean_eccentricity_interval() const {
  return mean_eccentricity_interval_;
}

inline Interval<Angle> OrbitalElements::mean_inclination_interval() const {
  return mean_inclination_interval_;
}

inline Interval<Angle>
OrbitalElements::mean_longitude_of_ascending_node_interval() const {
  return mean_longitude_of_ascending_node_interval_;
}

inline Interval<Angle> OrbitalElements::mean_argument_of_periapsis_interval()
    const {
  return mean_argument_of_periapsis_interval_;
}

inline Interval<Length> OrbitalElements::mean_periapsis_distance_interval()
    const {
  return mean_periapsis_distance_interval_;
}

inline Interval<Length> OrbitalElements::mean_apoapsis_distance_interval()
    const {
  return mean_apoapsis_distance_interval_;
}

inline Interval<Length> OrbitalElements::radial_distance_interval() const {
  return radial_distance_interval_;
}

template<typename PrimaryCentred>
std::vector<Length> OrbitalElements::RadialDistances(
    Trajectory<PrimaryCentred> const& trajectory) {
  std::vector<Length> radial_distances;
  auto const& discrete =
      dynamic_cast<physics::DiscreteTrajectory<PrimaryCentred> const&>(
          trajectory);
  radial_distances.reserve(discrete.size());
  DegreesOfFreedom<PrimaryCentred> const primary_dof{PrimaryCentred::origin,
                                                     PrimaryCentred::unmoving};
  for (auto const& [time, degrees_of_freedom] : discrete) {
    radial_distances.push_back(
        (degrees_of_freedom.position() - primary_dof.position()).Norm());
  }
  return radial_distances;
}

inline std::vector<OrbitalElements::EquinoctialElements> const&
OrbitalElements::osculating_equinoctial_elements() const {
  return osculating_equinoctial_elements_;
}

inline std::vector<OrbitalElements::EquinoctialElements> const&
OrbitalElements::mean_equinoctial_elements() const {
  return mean_equinoctial_elements_;
}

inline absl::StatusOr<Time> OrbitalElements::SiderealPeriod(
    std::function<EquinoctialElements(Instant const&)> const&
        equinoctial_elements,
    Instant const& t_min,
    Instant const& t_max) {
  Time const Δt = t_max - t_min;
  Instant const t0 = t_min + Δt / 2;
  int z = 0;
  Product<Angle, Square<Time>> const ʃ_λt_dt = AutomaticClenshawCurtis(
      [&equinoctial_elements, &t0, &z](
          Instant const& t) -> Product<Angle, Time> {
        // TODO(egg): Consider computing only λ.
        ++z;
        return equinoctial_elements(t).λ * (t - t0);
      },
      t_min,
      t_max,
      /*max_relative_error=*/1e-6,
      /*max_points=*/std::nullopt);
  LOG(ERROR) << "Automatic Clenshaw Curtis for sidereal period: " << z
             << " evaluations";
  return 2 * π * Radian * Pow<3>(Δt) / (12 * ʃ_λt_dt);
}

inline absl::StatusOr<std::vector<OrbitalElements::EquinoctialElements>>
OrbitalElements::MeanEquinoctialElements(
    std::function<EquinoctialElements(Instant const&)> const&
        equinoctial_elements,
    Instant const& t_min,
    Instant const& t_max,
    Time const& period) {
  // REMOVE BEFORE FLIGHT REWRITE THIS COMMENT
  // This function averages the elements in |osculating| over |period|.
  // For each |EquinoctialElements osculating_elements = osculating[i]| in
  // |osculating| such that |osculating_elements.t <= osculating.back().t|, let
  // |tᵢ = osculating_elements.t|; the result contains an
  // |EquinoctialElements mean_elements| with
  // |mean_elements.t = tᵢ + period / 2|.
  // For all э in the set of equinoctial elements {a, h, k, λ, p, q, pʹ, qʹ},
  // |mean_elements.э| is the integral of the osculating э from |tᵢ| to
  // |tᵢ + period|, divided by |period|.

  // Instead of computing the integral from |tᵢ| to |tᵢ + period| directly, we
  // precompute the integrals from |t_min| to each of the |tᵢ|.
  // The integral from |tᵢ| to |tᵢ + period| is then computed as the integral
  // from |tᵢ| to the last |tⱼ₋₁| before |tᵢ + period| (obtained by subtracting
  // the integrals to |tᵢ| and to |tⱼ₋₁|), plus the remainder integral from
  // |tⱼ₋₁| to |tᵢ + period| (a partial trapezoid on [tⱼ₋₁, tⱼ]).

  // Precompute the integrals from |t_min| to each of the |tᵢ|.


  int z = 0;

  /*
  for (Instant t = t_min; t <= t_max; t += period / 100) {
    auto const [_, a, h, k, λ, p, q, pʹ, qʹ] = equinoctial_elements(t);
    logger.Append("evaluatedOsculating",
                  std::tuple{t, a, h, k, λ, p, q, pʹ, qʹ},
                  mathematica::ExpressInSIUnits);
  }
  */

  using ODE =
      ExplicitFirstOrderOrdinaryDifferentialEquation<Instant,
                                                     Product<Length, Time>,
                                                     Time,
                                                     Time,
                                                     Product<Angle, Time>,
                                                     Time,
                                                     Time,
                                                     Time,
                                                     Time>;
  using namespace quantities::si;
  InitialValueProblem<ODE> problem{
      .equation =
          {
              .compute_derivative =
                  [&equinoctial_elements, &z](
                      Instant const& t,
                      ODE::DependentVariables const& y,
                      ODE::DependentVariableDerivatives& yʹ) {
                    ++z;
                    auto const [_, a, h, k, λ, p, q, pʹ, qʹ] =
                        equinoctial_elements(t);
                    // This rather ugly way of writing the assignments matter
                    // for performance: a temporary like {{a}, ..., {qʹ}} would
                    // do lots of allocations.
                    std::get<0>(yʹ)[0] = a;
                    std::get<1>(yʹ)[0] = h;
                    std::get<2>(yʹ)[0] = k;
                    std::get<3>(yʹ)[0] = λ;
                    std::get<4>(yʹ)[0] = p;
                    std::get<5>(yʹ)[0] = q;
                    std::get<6>(yʹ)[0] = pʹ;
                    std::get<7>(yʹ)[0] = qʹ;
                    return absl::OkStatus();
                  },
          },
      .initial_state = ODE::State(t_min,
                                  {{0 * Metre * Second},
                                   {0 * Second},
                                   {0 * Second},
                                   {0 * Radian * Second},
                                   {0 * Second},
                                   {0 * Second},
                                   {0 * Second},
                                   {0 * Second}}),
  };

  struct IntegratedEquinoctialElements {
    Instant t;
    // The integrals are from t_min to t_max.
    Product<Length, Time> ʃ_a_dt;
    Time ʃ_h_dt;
    Time ʃ_k_dt;
    Product<Angle, Time> ʃ_λ_dt;
    Time ʃ_p_dt;
    Time ʃ_q_dt;
    Time ʃ_pʹ_dt;
    Time ʃ_qʹ_dt;
  };
  std::vector<IntegratedEquinoctialElements> integrals;
  auto const append_state = [&integrals](ODE::State const& state) {
    Instant const& t = state.s.value;
    auto const& [ʃ_a_dt,
                 ʃ_h_dt,
                 ʃ_k_dt,
                 ʃ_λ_dt,
                 ʃ_p_dt,
                 ʃ_q_dt,
                 ʃ_pʹ_dt,
                 ʃ_qʹ_dt] = state.y;
    integrals.push_back(
        IntegratedEquinoctialElements{.t = t,
                                      .ʃ_a_dt = ʃ_a_dt.front().value,
                                      .ʃ_h_dt = ʃ_h_dt.front().value,
                                      .ʃ_k_dt = ʃ_k_dt.front().value,
                                      .ʃ_λ_dt = ʃ_λ_dt.front().value,
                                      .ʃ_p_dt = ʃ_p_dt.front().value,
                                      .ʃ_q_dt = ʃ_q_dt.front().value,
                                      .ʃ_pʹ_dt = ʃ_pʹ_dt.front().value,
                                      .ʃ_qʹ_dt = ʃ_qʹ_dt.front().value});
  };
  append_state(problem.initial_state);
  auto const instance =
      ExplicitLinearMultistepIntegrator<AdamsBashforthOrder6, ODE>()
          .NewInstance(problem,
                       append_state,
                       /*step=*/period / 24);
  RETURN_IF_ERROR(instance->Solve(t_max));
  LOG(ERROR) << z << " evaluations by integrator producing " << integrals.size()
             << " points";
  for (auto const& i : integrals) {
    auto const [_, a, h, k, λ, p, q, pʹ, qʹ] = equinoctial_elements(i.t);
  }

  auto const evaluate_integrals =
      [&integrals](Instant const& t) -> IntegratedEquinoctialElements {
    CHECK_LE(t, integrals.back().t);
    auto it = std::partition_point(
        integrals.begin(),
        integrals.end(),
        [&t](IntegratedEquinoctialElements const& elements) {
          return elements.t < t;
        });
    CHECK(it != integrals.end());
    IntegratedEquinoctialElements const& high = *it;
    if (it == integrals.begin()) {
      return high;
    } else {
      IntegratedEquinoctialElements const& low = *--it;
      double const α = (t - low.t) / (high.t - low.t);
      auto const interpolate = [α, &low, &high](auto const element) {
        return low.*element + α * (high.*element - low.*element);
      };
      return {.t = t,
              .ʃ_a_dt = interpolate(&IntegratedEquinoctialElements::ʃ_a_dt),
              .ʃ_h_dt = interpolate(&IntegratedEquinoctialElements::ʃ_h_dt),
              .ʃ_k_dt = interpolate(&IntegratedEquinoctialElements::ʃ_k_dt),
              .ʃ_λ_dt = interpolate(&IntegratedEquinoctialElements::ʃ_λ_dt),
              .ʃ_p_dt = interpolate(&IntegratedEquinoctialElements::ʃ_p_dt),
              .ʃ_q_dt = interpolate(&IntegratedEquinoctialElements::ʃ_q_dt),
              .ʃ_pʹ_dt = interpolate(&IntegratedEquinoctialElements::ʃ_pʹ_dt),
              .ʃ_qʹ_dt = interpolate(&IntegratedEquinoctialElements::ʃ_qʹ_dt)};
    }
  };

  // Now compute the averages.
  std::vector<EquinoctialElements> mean_elements;
  mean_elements.reserve(integrals.size());
  int j = 0;
  int y = 0;
  for (auto const& low : integrals) {
    RETURN_IF_STOPPED;
    if (low.t + period > integrals.back().t) {
      break;
    }
    auto const high = evaluate_integrals(low.t + period);
    mean_elements.emplace_back();
    mean_elements.back().t = low.t + period / 2;
    mean_elements.back().a = (high.ʃ_a_dt - low.ʃ_a_dt) / period;
    mean_elements.back().h = (high.ʃ_h_dt - low.ʃ_h_dt) / period;
    mean_elements.back().k = (high.ʃ_k_dt - low.ʃ_k_dt) / period;
    mean_elements.back().λ = (high.ʃ_λ_dt - low.ʃ_λ_dt) / period;
    mean_elements.back().p = (high.ʃ_p_dt - low.ʃ_p_dt) / period;
    mean_elements.back().q = (high.ʃ_q_dt - low.ʃ_q_dt) / period;
    mean_elements.back().pʹ = (high.ʃ_pʹ_dt - low.ʃ_pʹ_dt) / period;
    mean_elements.back().qʹ = (high.ʃ_qʹ_dt - low.ʃ_qʹ_dt) / period;
  }
  return mean_elements;
}

inline absl::StatusOr<std::vector<OrbitalElements::ClassicalElements>>
OrbitalElements::ToClassicalElements(
    std::vector<EquinoctialElements> const& equinoctial_elements) {
  std::vector<ClassicalElements> classical_elements;
  classical_elements.reserve(equinoctial_elements.size());
  for (auto const& equinoctial : equinoctial_elements) {
    RETURN_IF_STOPPED;
    double const tg_½i = Sqrt(Pow<2>(equinoctial.p) + Pow<2>(equinoctial.q));
    double const cotg_½i =
        Sqrt(Pow<2>(equinoctial.pʹ) + Pow<2>(equinoctial.qʹ));
    Angle const i =
        cotg_½i > tg_½i ? 2 * ArcTan(tg_½i) : 2 * ArcTan(1 / cotg_½i);
    Angle const Ω = cotg_½i > tg_½i ? ArcTan(equinoctial.p, equinoctial.q)
                                    : ArcTan(equinoctial.pʹ, equinoctial.qʹ);
    double const e = Sqrt(Pow<2>(equinoctial.h) + Pow<2>(equinoctial.k));
    Angle const ϖ = ArcTan(equinoctial.h, equinoctial.k);
    Angle const ω = ϖ - Ω;
    Angle const M = equinoctial.λ - ϖ;
    classical_elements.push_back(
        {.time = equinoctial.t,
         .semimajor_axis = equinoctial.a,
         .eccentricity = e,
         .inclination = i,
         .longitude_of_ascending_node = classical_elements.empty()
             ? Mod(Ω, 2 * π * Radian)
             : UnwindFrom(classical_elements.back().longitude_of_ascending_node,
                          Ω),
         .argument_of_periapsis = classical_elements.empty()
             ? Mod(ω, 2 * π * Radian)
             : UnwindFrom(classical_elements.back().argument_of_periapsis, ω),
         .mean_anomaly = classical_elements.empty()
             ? Mod(M, 2 * π * Radian)
             : UnwindFrom(classical_elements.back().mean_anomaly, M),
         .periapsis_distance = (1 - e) * equinoctial.a,
         .apoapsis_distance = (1 + e) * equinoctial.a});
  }
  return classical_elements;
}

inline absl::Status OrbitalElements::ComputePeriodsAndPrecession() {
  Time const Δt = mean_classical_elements_.back().time -
                  mean_classical_elements_.front().time;
  auto const Δt³ = Pow<3>(Δt);
  // We compute the mean rate (slope) of the mean anomaly M(t), the mean
  // argument of latitude u(t), and the longitude of the ascending node Ω(t).
  // On an interval [t̄ - Δt/2, t̄ + Δt/2], the slope of э is computed as
  //   ∫ (э(t) - э̄) (t - t̄) dt / ∫ (t - t̄)² dt;
  // this is the continuous analogue of a simple linear regression.
  // With ∫ (t - t̄)² dt = Δt³ / 12 and ∫ э̄ (t - t̄) = э̄ ∫ (t - t̄) dt = 0,
  // this simplifies to
  //   12 ∫ э(t) (t - t̄) dt / Δt³.
  // We first compute ∫ э(t) (t - t̄) dt for the three elements of interest.

  int z = 0;
  auto const interpolate_function_of_mean_classical_element =
      [this, &z](auto const f, Instant const& t) {
        CHECK_LE(t, mean_classical_elements_.back().time);
        auto it = std::partition_point(mean_classical_elements_.begin(),
                                       mean_classical_elements_.end(),
                                       [&t](ClassicalElements const& elements) {
                                         return elements.time < t;
                                       });
        ++z;
        ClassicalElements const& high = *it;
        if (it == mean_classical_elements_.begin()) {
          return f(high);
        } else {
          ClassicalElements const& low = *--it;
          double const α = (t - low.time) / (high.time - low.time);
          return f(low) + α * (f(high) - f(low));
        }
      };

  Instant const t̄ = mean_classical_elements_.front().time + Δt / 2;

  Product<Angle, Square<Time>> const ʃ_Mt_dt = AutomaticClenshawCurtis(
      [&interpolate_function_of_mean_classical_element, &t̄](Instant const& t) {
        return interpolate_function_of_mean_classical_element(
            [&t, &t̄](ClassicalElements const& elements) {
              return elements.mean_anomaly * (t - t̄);
            },
            t);
      },
      mean_classical_elements_.front().time,
      mean_classical_elements_.back().time,
      /*max_relative_error=*/1e-6,
      /*max_points=*/mean_classical_elements_.size());
  LOG(ERROR) << z << " for M";z=0;
  Product<Angle, Square<Time>> const ʃ_ut_dt = AutomaticClenshawCurtis(
      [&interpolate_function_of_mean_classical_element, &t̄](Instant const& t) {
        return interpolate_function_of_mean_classical_element(
            [&t, &t̄](ClassicalElements const& elements) {
              return (elements.argument_of_periapsis + elements.mean_anomaly) *
                     (t - t̄);
            },
            t);
      },
      mean_classical_elements_.front().time,
      mean_classical_elements_.back().time,
      /*max_relative_error=*/1e-6,
      /*max_points=*/mean_classical_elements_.size());
  LOG(ERROR) << z << " for u";z = 0;
  Product<Angle, Square<Time>> const ʃ_Ωt_dt = AutomaticClenshawCurtis(
      [&interpolate_function_of_mean_classical_element, &t̄](Instant const& t) {
        return interpolate_function_of_mean_classical_element(
            [&t, &t̄](ClassicalElements const& elements) {
              return elements.longitude_of_ascending_node * (t - t̄);
            },
            t);
      },
      mean_classical_elements_.front().time,
      mean_classical_elements_.back().time,
      /*max_relative_error=*/1e-6,
      /*max_points=*/mean_classical_elements_.size());
  LOG(ERROR) << z << " for Ω";

  // The periods are 2π over the mean rate of the relevant element; the nodal
  // precession is the mean rate of Ω.

  anomalistic_period_ = 2 * π * Radian * Δt³ / (12 * ʃ_Mt_dt);
  nodal_period_ = 2 * π * Radian * Δt³ / (12 * ʃ_ut_dt);
  nodal_precession_ = 12 * ʃ_Ωt_dt / Δt³;
  return absl::OkStatus();
}

inline absl::Status OrbitalElements::ComputeIntervals() {
  for (auto const& r : radial_distances_) {
    RETURN_IF_STOPPED;
    radial_distance_interval_.Include(r);
  }
  for (auto const& elements : mean_classical_elements_) {
    RETURN_IF_STOPPED;
    mean_semimajor_axis_interval_.Include(elements.semimajor_axis);
    mean_eccentricity_interval_.Include(elements.eccentricity);
    mean_inclination_interval_.Include(elements.inclination);
    mean_longitude_of_ascending_node_interval_.Include(
        elements.longitude_of_ascending_node);
    mean_argument_of_periapsis_interval_.Include(
        elements.argument_of_periapsis);
    mean_periapsis_distance_interval_.Include(elements.periapsis_distance);
    mean_apoapsis_distance_interval_.Include(elements.apoapsis_distance);
  }
  return absl::OkStatus();
}


}  // namespace internal_orbital_elements
}  // namespace astronomy
}  // namespace principia
