#pragma once

#include "astronomy/orbit_analyser.hpp"

#include "astronomy/orbit_recurrence.hpp"
#include "astronomy/orbital_elements.hpp"
#include "base/file.hpp"
#include "base/mod.hpp"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/apsides.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/kepler_orbit.hpp"
#include "quantities/astronomy.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbit_analyser {

using base::mod;
using geometry::AngularVelocity;
using geometry::Bivector;
using geometry::Displacement;
using geometry::OrientedAngleBetween;
using geometry::Position;
using geometry::Rotation;
using geometry::Sign;
using geometry::Vector;
using geometry::Velocity;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using physics::Body;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::ComputeApsides;
using physics::ComputeNodes;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::MasslessBody;
using physics::RigidMotion;
using physics::RigidTransformation;
using quantities::Abs;
using quantities::Angle;
using quantities::ArcTan;
using quantities::Cos;
using quantities::Infinity;
using quantities::Length;
using quantities::Pow;
using quantities::Product;
using quantities::Sin;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Square;
using quantities::Tan;
using quantities::Time;
using quantities::Variation;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

template<typename T>
T Average(std::vector<T> const& values) {
  T result;
  for (auto const& value : values) {
    result += (value - T{}) / values.size();
  }
  return result;
}

// Returns the element of {α + 2nπ| n ∈ ℤ} which is closest to |previous_angle|.
inline Angle UnwindFrom(Angle const& previous_angle, Angle const& α) {
  return α + std::nearbyint((previous_angle - α) / (2 * π * Radian)) * 2 * π *
                 Radian;
}

// |angles| should be sampled from a slowly-varying continuous function
// f: ℝ →  𝑆¹ = ℝ / 2π ℝ (specifically, consecutive angles should  differ by
// less than π).  Returns the corresponding sampling of the continuous g: ℝ → ℝ
// such that f = g mod 2π and f(0) = g(0).
inline std::vector<Angle> Unwind(std::vector<Angle> const& angles) {
  if (angles.empty()) {
    return angles;
  }
  std::vector<Angle> unwound_angles;
  unwound_angles.reserve(angles.size());
  unwound_angles.push_back(angles.front());
  for (int i = 1; i < angles.size(); ++i) {
    unwound_angles.push_back(UnwindFrom(unwound_angles.back(), angles[i]));
  }
  return unwound_angles;
}

template<typename T>
Difference<T> Variability(std::vector<T> const& values,
                          T const& nominal_value) {
  std::vector<Difference<T>> deviations;
  deviations.reserve(values.size());
  for (auto const& value : values) {
    deviations.push_back(Abs(value - nominal_value));
  }
  int n = 95 * values.size() / 100;
  std::nth_element(
      deviations.begin(), deviations.begin() + n, deviations.end());
  return deviations[n];
}

struct EquinoctialElements {
  Instant t;
  // See Broucke and Cefola (1972), On the equinoctial orbit elements.
  Length a;
  double h;
  double k;
  Angle λ;
  double p;
  double q;
  // pʹ and qʹ use the cotangent of the half-inclination instead of its tangent;
  // they are better suited to retrograde orbits.
  double pʹ;
  double qʹ;
};

std::vector<std::vector<double>> ElementsForLogging(
    std::vector<EquinoctialElements> const& elements_series) {
  std::vector<std::vector<double>> result;
  for (auto const& elements : elements_series) {
    result.push_back({(elements.t - J2000) / Second,
                      elements.a / Metre,
                      elements.h,
                      elements.k,
                      elements.λ / Radian,
                      elements.p,
                      elements.q,
                      elements.pʹ,
                      elements.qʹ});
  }
  return result;
}

template<typename PrimaryCentred>
std::vector<EquinoctialElements> OsculatingEquinoctialElements(
    DiscreteTrajectory<PrimaryCentred> const& trajectory,
    MassiveBody const& primary,
    Body const& secondary) {
  LOG(ERROR) << trajectory.Size();
  DegreesOfFreedom<PrimaryCentred> const primary_dof{
      PrimaryCentred::origin, Velocity<PrimaryCentred>{}};
  std::vector<EquinoctialElements> result;
  for (auto it = trajectory.Begin(); it != trajectory.End(); ++it) {
    auto const osculating_elements =
        KeplerOrbit<PrimaryCentred>(primary,
                                    secondary,
                                    it.degrees_of_freedom() - primary_dof,
                                    it.time())
            .elements_at_epoch();
    double const& e = *osculating_elements.eccentricity;
    Angle const& ϖ = *osculating_elements.longitude_of_periapsis;
    Angle const& Ω = osculating_elements.longitude_of_ascending_node;
    Angle const& M = *osculating_elements.mean_anomaly;
    Angle const& i = osculating_elements.inclination;
    double const tg_½i = Tan(i / 2);
    double const cotg_½i = 1 / tg_½i;
    result.push_back(
        {/*.t = */ it.time(),
         /*.a = */ *osculating_elements.semimajor_axis,
         /*.h = */ e * Sin(ϖ),
         /*.k = */ e * Cos(ϖ),
         /*.λ = */ result.empty() ? ϖ + M : UnwindFrom(result.back().λ, ϖ + M),
         /*.p = */ tg_½i * Sin(Ω),
         /*.q = */ tg_½i * Cos(Ω),
         /*.pʹ = */ cotg_½i * Sin(Ω),
         /*.qʹ = */ cotg_½i * Cos(Ω)});
  }
  LOG(ERROR) << result.size();
  return result;
}

// |equinoctial_elements| must contain at least 2 elements.
inline Time SiderealPeriod(
    std::vector<EquinoctialElements> equinoctial_elements) {
  LOG(ERROR) << equinoctial_elements.size();
  Time const Δt =
      equinoctial_elements.back().t - equinoctial_elements.front().t;
  Instant const t0 = equinoctial_elements.front().t + Δt / 2;
  Product<Angle, Square<Time>> ſ_λt_dt;

  for (auto previous = equinoctial_elements.begin(),
            it = equinoctial_elements.begin() + 1;
       it != equinoctial_elements.end();
       previous = it, ++it) {
    ſ_λt_dt += (it->λ * (it->t - t0) + previous->λ * (previous->t - t0)) / 2 *
               (it->t - previous->t);
  }
  return 2 * π * Radian * Pow<3>(Δt) / (12 * ſ_λt_dt);
}

// |osculating| must contain at least 2 elements.
// The resulting elements are averaged over one period, centred on t.
std::vector<EquinoctialElements> MeanEquinoctialElements(
    std::vector<EquinoctialElements> const& osculating,
    Time const& period) {
  Instant const& t_min = osculating.front().t;
  struct IntegratedEquinoctialElements {
    Instant t_max;
    // The integrals are from t_min to t_max.
    Product<Length, Time> ſ_a_dt;
    Time ſ_h_dt;
    Time ſ_k_dt;
    Product<Angle, Time> ſ_λ_dt;
    Time ſ_p_dt;
    Time ſ_q_dt;
    Time ſ_pʹ_dt;
    Time ſ_qʹ_dt;
  };
  std::vector<IntegratedEquinoctialElements> integrals;
  integrals.push_back({t_min});
  for (auto previous = osculating.begin(), it = osculating.begin() + 1;
       it != osculating.end();
       previous = it, ++it) {
    integrals.push_back(integrals.back());
    integrals.back().t_max = it->t;
    Time const dt = it->t - previous->t;
    integrals.back().ſ_a_dt += (it->a + previous->a) / 2 * dt;
    integrals.back().ſ_h_dt += (it->h + previous->h) / 2 * dt;
    integrals.back().ſ_k_dt += (it->k + previous->k) / 2 * dt;
    integrals.back().ſ_λ_dt += (it->λ + previous->λ) / 2 * dt;
    integrals.back().ſ_p_dt += (it->p + previous->p) / 2 * dt;
    integrals.back().ſ_q_dt += (it->q + previous->q) / 2 * dt;
    integrals.back().ſ_pʹ_dt += (it->pʹ + previous->pʹ) / 2 * dt;
    integrals.back().ſ_qʹ_dt += (it->qʹ + previous->qʹ) / 2 * dt;
  }
  std::vector<EquinoctialElements> result;
  int i = 0;
  for (auto first = integrals.begin(); first != integrals.end(); ++first) {
    Instant const t_first = first->t_max;
    while (integrals[i].t_max - t_first < period) {
      if (++i == integrals.size()) {
        return result;
      }
    }
    auto const& last = integrals[i - 1];
    // |element| should be a pointer to a member of |EquinoctialElements|;
    // Integrates from |last.t_max| to |t_first + period|.
    auto ſ = [i, &period, &t_first, &osculating](auto element) {
      auto const& next_osculating = osculating[i];
      auto const& last_osculating = osculating[i - 1];
      Time const Δt = next_osculating.t - last_osculating.t;
      Time const dt = t_first + period - last_osculating.t;
      auto const element_at_end =
          last_osculating.*element +
          (next_osculating.*element - last_osculating.*element) * (dt / Δt);
      return (last_osculating.*element + element_at_end) / 2 * dt;
    };
    result.emplace_back();
    result.back().t = t_first + period / 2;
    result.back().a =
        (last.ſ_a_dt - first->ſ_a_dt + ſ(&EquinoctialElements::a)) / period;
    result.back().h =
        (last.ſ_h_dt - first->ſ_h_dt + ſ(&EquinoctialElements::h)) / period;
    result.back().k =
        (last.ſ_k_dt - first->ſ_k_dt + ſ(&EquinoctialElements::k)) / period;
    result.back().λ =
        (last.ſ_λ_dt - first->ſ_λ_dt + ſ(&EquinoctialElements::λ)) / period;
    result.back().p =
        (last.ſ_p_dt - first->ſ_p_dt + ſ(&EquinoctialElements::p)) / period;
    result.back().q =
        (last.ſ_q_dt - first->ſ_q_dt + ſ(&EquinoctialElements::q)) / period;
    result.back().pʹ =
        (last.ſ_pʹ_dt - first->ſ_pʹ_dt + ſ(&EquinoctialElements::pʹ)) / period;
    result.back().qʹ =
        (last.ſ_qʹ_dt - first->ſ_qʹ_dt + ſ(&EquinoctialElements::qʹ)) / period;
  }
  return result;
}

// The classical Keplerian elements.
struct ClassicalElements {
  Instant t;
  Length a;
  double e;
  Angle Ω;
  Angle i;
  Angle ω;
  Angle M;
};

std::vector<ClassicalElements> ToClassicalElements(
    std::vector<EquinoctialElements> const& equinoctial_elements) {
  std::vector<ClassicalElements> result;
  for (auto const& equinoctial : equinoctial_elements) {
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
    result.push_back({equinoctial.t,
                      equinoctial.a,
                      e,
                      result.empty() ? Ω : UnwindFrom(result.back().Ω, Ω),
                      i,
                      result.empty() ? ω : UnwindFrom(result.back().ω, ω),
                      result.empty() ? M : UnwindFrom(result.back().M, M)});
  }
  return result;
}

// |elements| must contain at least 2 elements.
inline Time AnomalisticPeriod(std::vector<ClassicalElements> elements) {
  Time const Δt = elements.back().t - elements.front().t;
  Instant const t0 = elements.front().t + Δt / 2;
  Product<Angle, Square<Time>> ſ_Mt_dt;

  for (auto previous = elements.begin(), it = elements.begin() + 1;
       it != elements.end();
       previous = it, ++it) {
    ſ_Mt_dt += (it->M * (it->t - t0) + previous->M * (previous->t - t0)) / 2 *
               (it->t - previous->t);
  }
  return 2 * π * Radian * Pow<3>(Δt) / (12 * ſ_Mt_dt);
}

// |elements| must contain at least 2 elements.
inline Time NodalPeriod(std::vector<ClassicalElements> elements) {
  Time const Δt = elements.back().t - elements.front().t;
  Instant const t0 = elements.front().t + Δt / 2;
  Product<Angle, Square<Time>> ſ_ut_dt;

  auto previous = elements.begin();
  Angle previous_u = previous->ω + previous->M;
  for (auto it = elements.begin() + 1; it != elements.end();
       previous = it, ++it) {
    Angle const u = UnwindFrom(previous_u, it->ω + it->M);
    ſ_ut_dt += (u * (it->t - t0) + previous_u * (previous->t - t0)) / 2 *
               (it->t - previous->t);
    previous_u = u;
  }
  return 2 * π * Radian * Pow<3>(Δt) / (12 * ſ_ut_dt);
}

// |elements| must contain at least 2 elements.
inline Variation<Angle> NodalPrecession(
    std::vector<ClassicalElements> elements) {
  Time const Δt = elements.back().t - elements.front().t;
  Instant const t0 = elements.front().t + Δt / 2;
  Product<Angle, Square<Time>> ſ_Ωt_dt;

  for (auto previous = elements.begin(), it = elements.begin() + 1;
       it != elements.end();
       previous = it, ++it) {
    ſ_Ωt_dt += (it->Ω * (it->t - t0) + previous->Ω * (previous->t - t0)) / 2 *
               (it->t - previous->t);
  }
  return 12 * ſ_Ωt_dt / Pow<3>(Δt);
}

// |elements| must contain at least 2 elements.
inline Variation<Angle> ApsidalPrecession(
    std::vector<ClassicalElements> elements) {
  Time const Δt = elements.back().t - elements.front().t;
  Instant const t0 = elements.front().t + Δt / 2;
  Product<Angle, Square<Time>> ſ_ωt_dt;

  for (auto previous = elements.begin(), it = elements.begin() + 1;
       it != elements.end();
       previous = it, ++it) {
    ſ_ωt_dt += (it->ω * (it->t - t0) + previous->ω * (previous->t - t0)) / 2 *
               (it->t - previous->t);
  }
  return 12 * ſ_ωt_dt / Pow<3>(Δt);
}

template<typename Frame>
OrbitAnalyser<Frame>::OrbitAnalyser(
    not_null<Ephemeris<Frame>*> const ephemeris,
    not_null<RotatingBody<Frame> const*> const primary,
    not_null<RotatingBody<Frame> const*> sun,
    Instant const initial_time,
    DegreesOfFreedom<Frame> const initial_state,
    std::string name)
    : ephemeris_(ephemeris),
      primary_(primary),
      sun_(sun),
      name_(std::move(name)) {
  trajectory_.Append(initial_time, initial_state);
}

template<typename Frame>
OrbitAnalyser<Frame>::OrbitAnalyser(
    not_null<Ephemeris<Frame>*> const ephemeris,
    not_null<RotatingBody<Frame> const*> const primary,
    not_null<RotatingBody<Frame> const*> sun,
    DiscreteTrajectory<Frame> const& trajectory,
    std::string name)
    : ephemeris_(ephemeris),
      primary_(primary),
      sun_(sun),
      name_(std::move(name)) {
  for (auto it = trajectory.Begin(); it != trajectory.End(); ++it) {
    trajectory_.Append(it.time(), it.degrees_of_freedom());
  }

  Instant const t0 = trajectory.Begin().time();
  Time const osculating_year =
      *KeplerOrbit<Frame>(
           *sun_,
           *primary_,
           ephemeris_->trajectory(primary_)->EvaluateDegreesOfFreedom(t0) -
               ephemeris_->trajectory(sun_)->EvaluateDegreesOfFreedom(t0),
           t0)
           .elements_at_epoch()
           .period;
  // TODO(egg): We should probably have a way to compute the apsides and nodes
  // of a continuous trajectory.
  DiscreteTrajectory<Frame> primary_trajectory;
  for (double x = 0; x < 10; x += 1.0 / 1024) {
    Instant const t = t0 + x * osculating_year;
    ephemeris_->Prolong(t);
    primary_trajectory.Append(
        t, ephemeris_->trajectory(primary_)->EvaluateDegreesOfFreedom(t));
  }

  // We use the procedure from Lee (1995) as described in Allison and McEwen
  // (2000).
  // REMOVE BEFORE FLIGHT: cite the titles.
  // REMOVE BEFORE FLIGHT: do this with mean elements.
  DiscreteTrajectory<Frame> ascending_nodes;
  DiscreteTrajectory<Frame> descending_nodes;
  ComputeNodes(primary_trajectory.Begin(),
               primary_trajectory.End(),
               primary_->polar_axis(),
               ascending_nodes,
               descending_nodes);
  std::vector<Time> times_between_northward_equinoxes;
  {
    auto it = ascending_nodes.Begin();
    Instant previous_equinox = it.time();
    for (++it; it != ascending_nodes.End(); ++it) {
      times_between_northward_equinoxes.push_back(it.time() - previous_equinox);
      previous_equinox = it.time();
    }
  }
  tropical_year_ = Average(times_between_northward_equinoxes);

  DiscreteTrajectory<Frame> primary_aphelia;
  DiscreteTrajectory<Frame> primary_perihelia;
  ComputeApsides(*ephemeris_->trajectory(sun_),
                 primary_trajectory.Begin(),
                 primary_trajectory.End(),
                 primary_aphelia,
                 primary_perihelia);

  reference_perihelion_time_ = primary_perihelia.Begin().time();
  std::vector<Angle> adjusted_longitudes_of_perihelia;
  // REMOVE BEFORE FLIGHT: The reference direction here has to be node.
  for (auto it = primary_perihelia.Begin(); it != primary_perihelia.End();
       ++it) {
    Angle const argument_of_perihelion = OrientedAngleBetween(
        Vector<double, Frame>({1, 0, 0}),
        ephemeris_->trajectory(sun_)->EvaluatePosition(it.time()) -
            it.degrees_of_freedom().position(),
        Bivector<double, Frame>({0, 0, 1}));
    adjusted_longitudes_of_perihelia.push_back(
        argument_of_perihelion - 2 * π * Radian *
                                     (it.time() - reference_perihelion_time_) /
                                     tropical_year_);
  }
  longitude_of_perihelion_ = Average(Unwind(adjusted_longitudes_of_perihelia));
}

template<typename Frame>
void OrbitAnalyser<Frame>::Analyse() {
  Instant const t0 = trajectory_.Begin().time();
  ephemeris_->Prolong(t0);
  MasslessBody secondary;
  KeplerOrbit<Frame> const initial_osculating_orbit(
      *primary_,
      secondary,
      ephemeris_->trajectory(primary_)->EvaluateDegreesOfFreedom(t0) -
          trajectory_.Begin().degrees_of_freedom(),
      t0);
  KeplerianElements<Frame> const initial_osculating_elements =
      initial_osculating_orbit.elements_at_epoch();
  LOG(ERROR) << initial_osculating_elements;
  if (*initial_osculating_elements.eccentricity >= 1) {
    LOG(ERROR) << "Hyperbolic initial osculating orbit";
    return;
  }
  Time const& initial_osculating_period = *initial_osculating_elements.period;

  // TODO(egg): feed the uncertainties back into the tolerance so we don't use
  // overly small tolerances when the uncertainty is large anyway.
  // TODO(egg): consider using a fixed-step method.
  Length const length_tolerance =
      *initial_osculating_elements.semimajor_axis * 1e-7;
  Speed const speed_tolerance =
      2 * π * length_tolerance / initial_osculating_period;

  Ephemeris<Frame>::AdaptiveStepParameters parameters(
      EmbeddedExplicitRungeKuttaNyströmIntegrator<
          DormandالمكاوىPrince1986RKN434FM,
          Position<Frame>>(),
      /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
      length_tolerance,
      speed_tolerance);

  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 7 * Day,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 1 * JulianYear,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();
}

template<typename Frame>
void OrbitAnalyser<Frame>::RecomputeProperties() {
  LOG(ERROR) << "Analysing orbit of " << name_ << " over "
             << (trajectory_.last().time() - trajectory_.Begin().time()) / Day
             << " d = "
             << (trajectory_.last().time() - trajectory_.Begin().time()) /
                    JulianYear
             << " a";

  enum class PrimaryCentredTag { frame_tag };
  // The reference plane for orbital analysis is the xy plane.
  using PrimaryCentred = geometry::Frame<PrimaryCentredTag,
                                         PrimaryCentredTag::frame_tag,
                                         /*frame_is_inertial=*/true /*lies!*/>;
  BodyCentredNonRotatingDynamicFrame<Frame, PrimaryCentred> primary_centred(
      ephemeris_, primary_);

  DiscreteTrajectory<PrimaryCentred> primary_centred_trajectory;
  for (auto it = trajectory_.Begin(); it != trajectory_.End(); ++it) {
    primary_centred_trajectory.Append(
        it.time(),
        primary_centred.ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
  }

  auto orbital_elements = OrbitalElements::ForTrajectory(
          primary_centred_trajectory, *primary_, MasslessBody{})
          .ValueOrDie();

  auto const osculating_equinoctial_elements = OsculatingEquinoctialElements(
      primary_centred_trajectory, *primary_, MasslessBody{});
  sidereal_period_ = SiderealPeriod(osculating_equinoctial_elements);
  auto const osculating_classical_elements =
      ToClassicalElements(osculating_equinoctial_elements);
  auto const mean_equinoctial_elements = MeanEquinoctialElements(
      osculating_equinoctial_elements, sidereal_period_);
  auto const mean_classical_elements =
      ToClassicalElements(mean_equinoctial_elements);
  anomalistic_period_ = AnomalisticPeriod(mean_classical_elements);
  nodal_period_ = NodalPeriod(mean_classical_elements);
  nodal_precession_ = NodalPrecession(mean_classical_elements);
  apsidal_precession_ = ApsidalPrecession(mean_classical_elements);

  {
    base::OFStream file(SOLUTION_DIR / (name_ + "_elements"));
    file << mathematica::Assign(
        name_ + "osculatingEquinoctialElements",
        ElementsForLogging(osculating_equinoctial_elements));
    file << mathematica::Assign(name_ + "meanEquinoctialElements",
                                ElementsForLogging(mean_equinoctial_elements));
  }

  auto const recurrence = OrbitRecurrence::ClosestRecurrence(
      nodal_period_, nodal_precession_, *primary_, /*max_abs_Cᴛₒ=*/50);
  Angle const Δλᴇ = recurrence.equatorial_shift();
  Angle const δʀ = recurrence.base_interval();

  DiscreteTrajectory<PrimaryCentred> ascending_nodes;
  DiscreteTrajectory<PrimaryCentred> descending_nodes;
  ComputeNodes(primary_centred_trajectory.Begin(),
               primary_centred_trajectory.End(),
               Vector<double, PrimaryCentred>({0, 0, 1}),
               ascending_nodes,
               descending_nodes);

  LOG(ERROR) << ascending_nodes.Size() << " ascending nodes";

  std::vector<Angle> reduced_terrestrial_longitudes_of_ascending_nodes;
  std::vector<Angle> mean_solar_times_of_ascending_nodes;
  int k;
  int i = 0;
  for (auto node = ascending_nodes.Begin(); node != ascending_nodes.End();
       ++node, ++i) {
    Angle const Ω =
        (node.degrees_of_freedom().position() - PrimaryCentred::origin)
            .coordinates()
            .ToSpherical()
            .longitude;
    Instant const t = node.time();
    // TODO(egg): this assumes earthlike sign for the longitude.
    Angle const λ = Ω - (primary_->AngleAt(t) + π / 2 * Radian);
    // k is the initial number of equatorial shifts from the interval of
    // longitudes [0, δʀ].
    // REMOVE BEFORE FLIGHT: this is probably wrong for Δλᴇ > 0; further,
    // we probably want to reduce on the grid interval [0, δ].
    if (i == 0) {
      k = Sign(Δλᴇ) * std::floor(λ / Abs(Δλᴇ));
    }
    Angle const λ_reduced =
        quantities::Mod(λ - (i + k) * Δλᴇ + π * Radian, 2 * π * Radian) -
        π * Radian;
    reduced_terrestrial_longitudes_of_ascending_nodes.push_back(
        i == 0 ? λ_reduced
               : UnwindFrom(
                     reduced_terrestrial_longitudes_of_ascending_nodes.back(),
                     λ_reduced));

    // Ignoring the error bars on the mean sun, effectively making it
    // conventional.
    Angle const mean_solar_time =
        quantities::Mod(
            (node.degrees_of_freedom().position() - PrimaryCentred::origin)
                    .coordinates()
                    .ToSpherical()
                    .longitude -
                (longitude_of_perihelion_ +
                 2 * π * Radian * (node.time() - reference_perihelion_time_) /
                     tropical_year_) +
                π * Radian,
            2 * π * Radian) -
        π * Radian;
    mean_solar_times_of_ascending_nodes.push_back(mean_solar_time);
  }

  auto const τ = Unwind(mean_solar_times_of_ascending_nodes);
  Angle const mean_τ = Average(τ);

  auto const λ0 = Average(reduced_terrestrial_longitudes_of_ascending_nodes);

  LOG(ERROR) << "--- General parameters ---";
  LOG(ERROR) << u8"T* = " << sidereal_period_ / Second << " s";
  LOG(ERROR) << u8"Td = " << nodal_period_ / Second << " s";
  LOG(ERROR) << u8"T = " << anomalistic_period_ / Second << " s";
  LOG(ERROR) << u8"Ω′ = " << nodal_precession_ / (Degree / Day) << u8"°/d";
  LOG(ERROR) << u8"ω′ = " << apsidal_precession_ / (Degree / Day) << u8"°/d";
  LOG(ERROR) << "i = ";
  LOG(ERROR) << "OrbitalElements";
  LOG(ERROR) << orbital_elements.sidereal_period() / Second << " s";
  LOG(ERROR) << orbital_elements.nodal_period() / Second << " s";
  LOG(ERROR) << orbital_elements.anomalistic_period() / Second << " s";
  LOG(ERROR) << u8"Ω′ = "
             << orbital_elements.nodal_precession() / (Degree / Day) << u8"°/d";
  LOG(ERROR) << u8"a = ["
             << orbital_elements.mean_semimajor_axis().min / (Kilo(Metre))
             << u8" km,";
  LOG(ERROR) << u8"     "
             << orbital_elements.mean_semimajor_axis().max / (Kilo(Metre))
             << u8" km]";
  LOG(ERROR) << u8"e = [" << orbital_elements.mean_eccentricity().min << ",";
  LOG(ERROR) << u8"     " << orbital_elements.mean_eccentricity().max << "]";
  LOG(ERROR) << u8"i = [" << orbital_elements.mean_inclination().min / Degree
             << u8"°,";
  LOG(ERROR) << u8"     " << orbital_elements.mean_inclination().max / Degree
             << u8"°]";
  LOG(ERROR) << u8"Ω = ["
             << orbital_elements.mean_longitude_of_ascending_node().min / Degree
             << u8"°,";
  LOG(ERROR) << u8"     "
             << orbital_elements.mean_longitude_of_ascending_node().max / Degree
             << u8"°]";
  LOG(ERROR) << u8"ω = ["
             << orbital_elements.mean_argument_of_periapsis().min / Degree
             << u8"°,";
  LOG(ERROR) << u8"     "
             << orbital_elements.mean_argument_of_periapsis().max / Degree
             << u8"°]";
  LOG(ERROR) << "----";
  LOG(ERROR) << "--- Orbit with respect to the Earth ---";
  LOG(ERROR) << "- Phasing -";
  LOG(ERROR) << "N_To / C_To = " << recurrence.number_of_revolutions() << " / "
             << recurrence.Cᴛₒ();
  LOG(ERROR) << u8"[ν0 ; DTo ; CTo] = [" << recurrence.νₒ() << " ; "
             << recurrence.Dᴛₒ() << " ; " << recurrence.Cᴛₒ() << "]";
  LOG(ERROR) << u8"𝕃 = NTo Td = "
             << recurrence.number_of_revolutions() * nodal_period_ / Day
             << " d";

  LOG(ERROR) << u8"ΔλE = " << Δλᴇ / Degree << u8"°";
  LOG(ERROR) << u8"δ = " << recurrence.grid_interval() / Degree << u8"°";
  LOG(ERROR) << u8"ETo* = " << recurrence.subcycle();

  LOG(ERROR) << "- Ground track -";
  // TODOegg): extremal latitudes (& nominal inclination);
  LOG(ERROR) << u8"λ0 =" << λ0 / Degree << u8"°";
  LOG(ERROR) << u8"   ± "
             << Variability(reduced_terrestrial_longitudes_of_ascending_nodes,
                            λ0) /
                    Degree
             << u8"° (95%)";

  if (recurrence.νₒ() == 1 && recurrence.Cᴛₒ() == 1) {
    LOG(ERROR) << "- Geosynchronous orbit: central longitude -";
  }

  LOG(ERROR) << "- Orbit freezing -";

  LOG(ERROR) << "--- Orbit with respect to the Sun, heliosynchronicity ---";

  // REMOVE BEFORE FLIGHT: There probably should be a sign here to turn the
  // tropical year into Ω'S.
  auto const ΔtS =
      2 * π * Radian / (nodal_precession_ - (2 * π * Radian / tropical_year_));
  LOG(ERROR) << u8"ΔtS = " << ΔtS / Day << " d";

  LOG(ERROR) << u8"τNA = " << mean_τ / Degree << u8"° = "
             << 12 + (mean_τ * 24 / (2 * π * Radian)) << u8" h";
  LOG(ERROR) << u8"       ± "
             << Variability(τ, mean_τ) * (24 / (2 * π * Radian)) << " h (95 %)";
}

}  // namespace internal_orbit_analyser
}  // namespace astronomy
}  // namespace principia
