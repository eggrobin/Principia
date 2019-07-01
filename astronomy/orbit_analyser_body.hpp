#pragma once

#include "astronomy/orbit_analyser.hpp"
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
using quantities::AverageOfCorrelated;
using quantities::Cos;
using quantities::Infinity;
using quantities::Length;
using quantities::LinearRegression;
using quantities::Pow;
using quantities::Product;
using quantities::Sin;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Square;
using quantities::Tan;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

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
    result.push_back({(elements.t - elements_series.front().t) / Second,
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
                      Ω,
                      i,
                      ω,
                      result.empty() ? M : UnwindFrom(result.back().M, M)});
  }
  return result;
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
  tropical_year_ = AverageOfCorrelated(times_between_northward_equinoxes);

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
                                     tropical_year_.measured_value);
  }
  longitude_of_perihelion_ =
      AverageOfCorrelated(Unwind(adjusted_longitudes_of_perihelia));
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

  enum class PrimaryTag { normal, sideways };
  // The origin of the reference frame is the centre of mass of |*primary_|.
  // The axes are those of |Frame|.
  // The reference plane for orbital analysis is the xy plane.
  using PrimaryCentred = geometry::Frame<PrimaryTag,
                                         PrimaryTag::normal,
                                         /*frame_is_inertial=*/true /*lies!*/>;
  Vector<double, PrimaryCentred> x({1, 0, 0});
  Bivector<double, PrimaryCentred> z({0, 0, 1});

  // Same as PrimaryCentredNonRotating, with different axes:
  // - the x axis is the same as the one of PrimaryCentred;
  // - the y axis is the z axis of PrimaryCentred;
  // - the z axis is opposite the y axis of PrimaryCentred.
  using PrimaryCentredSideways = geometry::Frame<PrimaryTag,
                                                 PrimaryTag::sideways,
                                                 /*frame_is_inertial=*/false>;

  RigidMotion<PrimaryCentred, PrimaryCentredSideways> tip(
      RigidTransformation<PrimaryCentred, PrimaryCentredSideways>(
          PrimaryCentred::origin,
          PrimaryCentredSideways::origin,
          Rotation<PrimaryCentred, PrimaryCentredSideways>(x, z, x * z)
              .Forget()),
      /*angular_velocity_of_to_frame=*/AngularVelocity<PrimaryCentred>{},
      /*velocity_of_to_frame_origin=*/
      Velocity<PrimaryCentred>{});

  BodyCentredNonRotatingDynamicFrame<Frame, PrimaryCentred> primary_centred(
      ephemeris_, primary_);

  DiscreteTrajectory<PrimaryCentred> primary_centred_trajectory;
  DiscreteTrajectory<PrimaryCentredSideways>
      sideways_primary_centred_trajectory;
  for (auto it = trajectory_.Begin(); it != trajectory_.End(); ++it) {
    primary_centred_trajectory.Append(
        it.time(),
        primary_centred.ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
    sideways_primary_centred_trajectory.Append(
        it.time(), tip(primary_centred_trajectory.last().degrees_of_freedom()));
  }

  auto const osculating_equinoctial_elements = OsculatingEquinoctialElements(
      primary_centred_trajectory, *primary_, MasslessBody{});
  auto const sidereal_period = SiderealPeriod(osculating_equinoctial_elements);
  LOG(ERROR) << "sidereal period by integration = " << sidereal_period;
  auto const mean_equinoctial_elements =
      MeanEquinoctialElements(osculating_equinoctial_elements, sidereal_period);
  auto const mean_classical_elements =
      ToClassicalElements(mean_equinoctial_elements);

  {
    base::OFStream file(SOLUTION_DIR / (name_ + "_elements"));
    file << mathematica::Assign(
        name_ + "osculatingEquinoctialElements",
        ElementsForLogging(osculating_equinoctial_elements));
    file << mathematica::Assign(
        name_ + "meanEquinoctialElements",
        ElementsForLogging(mean_equinoctial_elements));
  }

  // REMOVE BEFORE FLIGHT: we should pick a reference direction orthogonal to
  // the orbital plane here, to avoid issues with orbits in the xz plane.

  // By computing the "nodes" with respect to the xz plane, i.e., the crossings
  // of the xz plane, we compute the points at which the projection of the orbit
  // onto the reference plane xy passes the fixed reference direction x.
  // While these points themselves are largely meaningless, the time between
  // them is the sidereal period.
  DiscreteTrajectory<PrimaryCentredSideways> xz_ascensions;
  DiscreteTrajectory<PrimaryCentredSideways> xz_descents;
  ComputeNodes(sideways_primary_centred_trajectory.Begin(),
               sideways_primary_centred_trajectory.End(),
               tip.orthogonal_map()(z * x),
               xz_ascensions,
               xz_descents);
  std::vector<Time> times_between_xz_ascensions;
  auto ascension = xz_ascensions.Begin();
  Instant previous_time = ascension.time();
  if (ascension != xz_ascensions.End()) {
    ++ascension;
  }
  for (; ascension != xz_ascensions.End();
       previous_time = ascension.time(), ++ascension) {
    times_between_xz_ascensions.push_back(ascension.time() - previous_time);
  }
  sidereal_period_ = AverageOfCorrelated(times_between_xz_ascensions);
  Time const sidereal_rotation_period =
      2 * π * Radian / primary_->angular_frequency();

  DiscreteTrajectory<PrimaryCentred> ascending_nodes;
  DiscreteTrajectory<PrimaryCentred> descending_nodes;
  ComputeNodes(primary_centred_trajectory.Begin(),
               primary_centred_trajectory.End(),
               Vector<double, PrimaryCentred>({0, 0, 1}),
               ascending_nodes,
               descending_nodes);

  LOG(ERROR) << ascending_nodes.Size() << " ascending nodes";

  std::vector<Instant> times_of_ascending_nodes;
  std::vector<Time> times_between_ascending_nodes;
  std::vector<Angle> longitudes_of_ascending_nodes;
  std::vector<Angle> inclinations_at_ascending_nodes;
  std::vector<Angle> terrestrial_longitudes_of_ascending_nodes;
  std::vector<Angle> true_solar_times_of_ascending_nodes;
  std::vector<Angle> mean_solar_times_of_ascending_nodes;
  for (auto node = ascending_nodes.Begin(); node != ascending_nodes.End();
       ++node) {
    // We do not construct |KeplerianElements|: we only need the longitude of
    // the ascending node, and we are at the ascending node so the computation
    // is trivial.
    Angle const Ω = OrientedAngleBetween(
        x, node.degrees_of_freedom().position() - PrimaryCentred::origin, z);
    if (!times_of_ascending_nodes.empty()) {
      times_between_ascending_nodes.push_back(node.time() -
                                              times_of_ascending_nodes.back());
    }
    times_of_ascending_nodes.push_back(node.time());
    longitudes_of_ascending_nodes.push_back(Ω);

    // REMOVE BEFORE FLIGHT: we should be able to orthogonalize vectors against
    // bivectors and vice-versa.
    // 0 for noon.
    Angle const true_solar_time = OrientedAngleBetween(
        (primary_centred
             .ToThisFrameAtTime(node.time())(
                 ephemeris_->trajectory(sun_)->EvaluateDegreesOfFreedom(
                     node.time()))
             .position() -
         PrimaryCentred::origin)
            .OrthogonalizationAgainst(
                Vector<double, PrimaryCentred>({0, 0, 1})),
        node.degrees_of_freedom().position() - PrimaryCentred::origin,
        z);
    true_solar_times_of_ascending_nodes.push_back(true_solar_time);

    // Ignoring the error bars on the mean sun, effectively making it
    // conventional.
    Angle const mean_solar_time =
        quantities::Mod(
            (node.degrees_of_freedom().position() - PrimaryCentred::origin)
                    .coordinates()
                    .ToSpherical()
                    .longitude -
                (longitude_of_perihelion_.measured_value +
                 2 * π * Radian * (node.time() - reference_perihelion_time_) /
                     tropical_year_.measured_value) +
                π * Radian,
            2 * π * Radian) -
        π * Radian;
    mean_solar_times_of_ascending_nodes.push_back(mean_solar_time);

    inclinations_at_ascending_nodes.push_back(AngleBetween(
        z,
        Wedge(node.degrees_of_freedom().position() - PrimaryCentred::origin,
              node.degrees_of_freedom().velocity())));
  }
  nodal_precession_ = LinearRegression(times_of_ascending_nodes,
                                       Unwind(longitudes_of_ascending_nodes))
                          .slope;
  nodal_period_ = AverageOfCorrelated(times_between_ascending_nodes);

  auto const tsv = Unwind(true_solar_times_of_ascending_nodes);
  MeasurementResult<Angle> const mean_tsv = AverageOfCorrelated(tsv); /*
  LOG(ERROR) << u8"TSV_NA = " << mean_tsv / Degree << u8"° = "
             << 12 + (mean_tsv * 24 / (2 * π * Radian)) << u8" h";
  LOG(ERROR) << u8"       ± "
             << Variability(tsv, mean_tsv.measured_value) *
                    (24 / (2 * π * Radian))
             << " h (95 %)";*/

  auto const τ = Unwind(mean_solar_times_of_ascending_nodes);
  MeasurementResult<Angle> const mean_τ = AverageOfCorrelated(τ);

  // TODO(egg): Consider factoring this out.
  std::vector<Angle> inclinations_at_extremal_latitudes;
  std::vector<Angle> all_latitudes;
  std::vector<AngularFrequency> all_latitude_rates;
  std::vector<Instant> all_times;
  {
    auto const latitude = [](Position<PrimaryCentred> q) -> Angle {
      return (q - PrimaryCentred::origin).coordinates().ToSpherical().latitude;
    };
    auto const latitude_rate =
        [](DegreesOfFreedom<PrimaryCentred> dof) -> AngularFrequency {
      Displacement<PrimaryCentred> const r =
          dof.position() - PrimaryCentred::origin;
      Vector<double, PrimaryCentred> const celestial_north({0, 0, 1});
      Vector<double, PrimaryCentred> const local_north =
          Normalize(celestial_north.OrthogonalizationAgainst(r));
      return InnerProduct(dof.velocity(), local_north) * Radian / r.Norm();
    };
    auto it = primary_centred_trajectory.Begin();
    Instant previous_time = it.time();
    Angle previous_latitude = latitude(it.degrees_of_freedom().position());
    AngularFrequency previous_latitude_rate =
        latitude_rate(it.degrees_of_freedom());
    for (++it; it != primary_centred_trajectory.End(); ++it) {
      Angle const new_latitude = latitude(it.degrees_of_freedom().position());
      AngularFrequency const new_latitude_rate =
          latitude_rate(it.degrees_of_freedom());
      all_times.push_back(it.time());
      all_latitudes.push_back(new_latitude);
      all_latitude_rates.push_back(new_latitude_rate);
      if (geometry::Sign(new_latitude_rate) !=
          geometry::Sign(previous_latitude_rate)) {
        numerics::Hermite3<Instant, Angle> interpolated_latitude(
            {previous_time, it.time()},
            {previous_latitude, new_latitude},
            {previous_latitude_rate, new_latitude_rate});
        Instant extremum_time;
        int valid_extrema = 0;
        for (Instant const& extremum : interpolated_latitude.FindExtrema()) {
          if (extremum >= previous_time && extremum <= it.time()) {
            extremum_time = extremum;
            ++valid_extrema;
          }
        }
        if (valid_extrema != 1) {
          extremum_time = geometry::Barycentre<Instant, AngularFrequency>(
              {it.time(), previous_time},
              {previous_latitude_rate, -new_latitude_rate});
        }
        inclinations_at_extremal_latitudes.push_back(Abs(latitude(
            primary_centred_trajectory.EvaluatePosition(extremum_time))));
        if (geometry::Sign(
                InnerProduct(Wedge((primary_centred_trajectory.EvaluatePosition(
                                        extremum_time) -
                                    PrimaryCentred::origin),
                                   primary_centred_trajectory.EvaluateVelocity(
                                       extremum_time)),
                             z))
                .Negative()) {
          inclinations_at_extremal_latitudes.back() =
              π * Radian - inclinations_at_extremal_latitudes.back();
        }
      }
      previous_time = it.time();
      previous_latitude = new_latitude;
      previous_latitude_rate = new_latitude_rate;
    }
  }
  // TODO(egg): this would need special handling for retrograde orbits; more
  // worryingly it is unsound for polar orbits.
  inclination_ = AverageOfCorrelated(inclinations_at_extremal_latitudes);
  LOG(ERROR) << "i = "
             << AverageOfCorrelated(inclinations_at_ascending_nodes) / Degree
             << u8"° (i☊)";
  DiscreteTrajectory<Frame> periapsides;
  DiscreteTrajectory<Frame> apoapsides;
  ComputeApsides(*ephemeris_->trajectory(primary_),
                 trajectory_.Begin(),
                 trajectory_.End(),
                 apoapsides,
                 periapsides);

  LOG(ERROR) << periapsides.Size() << " apparent periapsides";
  LOG(ERROR) << apoapsides.Size() << " apoapsides";

  std::vector<Instant> times_of_periapsides;
  std::vector<Instant> times_of_apoapsides;
  std::vector<Time> times_between_periapsides;
  std::vector<Angle> arguments_of_periapsides;
  std::vector<Length> periapsis_distances;
  std::vector<Length> apoapsis_distances;

  for (std::optional<DiscreteTrajectory<Frame>::Iterator> previous_periapsis;
       ;) {
    // Eliminate spurious periapsides at low eccentricities by taking the lowest
    // periapsis over 1.75 nodal periods (which should be enough to cover 1
    // apsidal period, but not enough to cover 2). REMOVE BEFORE FLIGHT: Make
    // sure that this really cannot cover 2 apsidal periods.
    Length periapsis_distance = Infinity<Length>();
    auto tentative_periapsis = previous_periapsis.value_or(periapsides.Begin());
    if (previous_periapsis.has_value()) {
      ++tentative_periapsis;
    }
    auto periapsis = tentative_periapsis;
    for (; tentative_periapsis != periapsides.End() &&
           tentative_periapsis.time() -
                   previous_periapsis.value_or(periapsides.Begin()).time() <=
               nodal_period_.measured_value * 1.75;
         ++tentative_periapsis) {
      Length const tentative_periapsis_distance =
          (ephemeris_->trajectory(primary_)->EvaluatePosition(
               tentative_periapsis.time()) -
           tentative_periapsis.degrees_of_freedom().position())
              .Norm();
      if (tentative_periapsis_distance < periapsis_distance) {
        periapsis_distance = tentative_periapsis_distance;
        periapsis = tentative_periapsis;
      }
    }
    if (tentative_periapsis == periapsides.End()) {
      // We have not been able to look at a whole 1.75 nodal periods after the
      // last periapsis; ignore trailing periapsides, as they may be spurious.
      break;
    }

    // TODO(egg): We could probably do something a lot more efficient, because
    // we know that we are at the periapsis, and we only need the argument of
    // periapsis.
    auto const elements =
        KeplerOrbit<Frame>(
            *primary_,
            MasslessBody{},
            periapsis.degrees_of_freedom() -
                ephemeris_->trajectory(primary_)->EvaluateDegreesOfFreedom(
                    periapsis.time()),
            periapsis.time())
            .elements_at_epoch();
    Angle ω = quantities::Mod(
        *elements.argument_of_periapsis + *elements.true_anomaly,
        2 * π * Radian);
    if (!times_of_periapsides.empty()) {
      times_between_periapsides.push_back(periapsis.time() -
                                          times_of_periapsides.back());
    }
    times_of_periapsides.push_back(periapsis.time());
    arguments_of_periapsides.push_back(ω);
    periapsis_distances.push_back(periapsis_distance);

    previous_periapsis = periapsis;
  }
  LOG(ERROR) << periapsis_distances.size() << " true periapsides";

  for (std::optional<DiscreteTrajectory<Frame>::Iterator> previous_apoapsis;;) {
    // Same as above, for apoapsides.
    Length apoapsis_distance = Infinity<Length>();
    auto tentative_apoapsis = previous_apoapsis.value_or(apoapsides.Begin());
    if (previous_apoapsis.has_value()) {
      ++tentative_apoapsis;
    }
    auto apoapsis = tentative_apoapsis;
    for (; tentative_apoapsis != apoapsides.End() &&
           tentative_apoapsis.time() -
                   previous_apoapsis.value_or(apoapsides.Begin()).time() <=
               nodal_period_.measured_value * 1.75;
         ++tentative_apoapsis) {
      Length const tentative_apoapsis_distance =
          (ephemeris_->trajectory(primary_)->EvaluatePosition(
               tentative_apoapsis.time()) -
           tentative_apoapsis.degrees_of_freedom().position())
              .Norm();
      if (tentative_apoapsis_distance < apoapsis_distance) {
        apoapsis_distance = tentative_apoapsis_distance;
        apoapsis = tentative_apoapsis;
      }
    }
    if (tentative_apoapsis == apoapsides.End()) {
      // We have not been able to look at a whole 1.75 nodal periods after the
      // last apoapsis; ignore trailing apoapsides, as they may be spurious.
      break;
    }
    times_of_apoapsides.push_back(apoapsis.time());
    apoapsis_distances.push_back(apoapsis_distance);

    previous_apoapsis = apoapsis;
  }
  LOG(ERROR) << apoapsis_distances.size() << " true apoapsides";

  apsidal_precession_ =
      LinearRegression(times_of_periapsides, Unwind(arguments_of_periapsides))
          .slope;
  anomalistic_period_ = AverageOfCorrelated(times_between_periapsides);
  periapsis_distance_ = AverageOfCorrelated(periapsis_distances);
  apoapsis_distance_ = AverageOfCorrelated(apoapsis_distances);
  eccentricity_ = 1 - 2 / (apoapsis_distance_ / periapsis_distance_ + 1);
  auto const ω = AverageOfCorrelated(arguments_of_periapsides);

  quantities::Product<Length, Time> ſ_a_dt;
  quantities::Product<Length, Time> ſ_r_pe_dt;
  quantities::Product<Length, Time> ſ_r_ap_dt;
  Time ſ_e_dt;
  Time ſ_e_cos_ω_dt;
  Time ſ_e_sin_ω_dt;
  Time ſ_log_e_dt;
  Time ſ_e⁻¹_dt;
  quantities::Product<Angle, Time> ſ_i_dt;
  std::vector<Instant> times;
  std::vector<Angle> mean_anomalies;
  std::vector<Angle> mean_longitudes;
  std::vector<Angle> mean_arguments_of_latitude;
  {
    std::optional<Instant> previous_time;
    std::optional<Length> previous_a;
    std::optional<Length> previous_r_pe;
    std::optional<Length> previous_r_ap;
    std::optional<double> previous_e;
    std::optional<Angle> previous_ω;
    std::optional<Angle> previous_i;
    for (auto it = trajectory_.Begin(); it != trajectory_.End(); ++it) {
      auto const elements =
          KeplerOrbit<Frame>(
              *primary_,
              MasslessBody{},
              it.degrees_of_freedom() -
                  ephemeris_->trajectory(primary_)->EvaluateDegreesOfFreedom(
                      it.time()),
              it.time())
              .elements_at_epoch();
      times.push_back(it.time());
      mean_anomalies.push_back(*elements.mean_anomaly);
      mean_longitudes.push_back(*elements.longitude_of_periapsis +
                                *elements.mean_anomaly);
      mean_arguments_of_latitude.push_back(*elements.argument_of_periapsis +
                                           *elements.mean_anomaly);
      if (previous_time.has_value()) {
        Time const Δt = it.time() - *previous_time;
        ſ_a_dt += (*previous_a + *elements.semimajor_axis) / 2 * Δt;
        ſ_r_pe_dt += (*previous_r_pe + *elements.periapsis_distance) / 2 * Δt;
        ſ_r_ap_dt += (*previous_r_ap + *elements.apoapsis_distance) / 2 * Δt;
        ſ_e_cos_ω_dt +=
            (*previous_e * Cos(*previous_ω) +
             *elements.eccentricity * Cos(*elements.argument_of_periapsis)) /
            2 * Δt;
        ſ_e_sin_ω_dt += (*previous_e * quantities::Sin(*previous_ω) +
                         *elements.eccentricity *
                             quantities::Sin(*elements.argument_of_periapsis)) /
                        2 * Δt;
        ſ_e_dt += (*previous_e + *elements.eccentricity) / 2 * Δt;
        ſ_log_e_dt +=
            (std::log(*previous_e) + std::log(*elements.eccentricity)) / 2 * Δt;
        ſ_e⁻¹_dt += (1 / *previous_e + 1 / *elements.eccentricity) / 2 * Δt;
        // TODO(egg): We should probably unwind.
        ſ_i_dt += (*previous_i + elements.inclination) / 2 * Δt;
      }
      previous_time = it.time();
      previous_a = elements.semimajor_axis;
      previous_r_pe = elements.periapsis_distance;
      previous_r_ap = elements.apoapsis_distance;
      previous_e = elements.eccentricity;
      previous_ω = elements.argument_of_periapsis;
      previous_i = elements.inclination;
    }
  }

  // (7.41).
  MeasurementResult<double> const daily_recurrence_frequency =
      (2 * π * Radian / nodal_period_) /
      (primary_->angular_frequency() - nodal_precession_);

  // 11.7.2.
  double smallest_fraction = std::numeric_limits<double>::infinity();
  int cycle_days;
  int nto;
  int cto;
  for (int j = 1; j < 50; ++j) {
    MeasurementResult<double> κ_j = daily_recurrence_frequency * j;
    double const abs_κ_j = std::abs(κ_j.measured_value);
    double const fraction = std::abs(abs_κ_j - std::nearbyint(abs_κ_j));
    if (fraction < smallest_fraction) {
      cycle_days = j;
      smallest_fraction = fraction;
      nto = std::nearbyint(abs_κ_j);
      cto = j;
    }
  }

  int const ν0 = std::nearbyint(daily_recurrence_frequency.measured_value);
  int const dto = nto - ν0 * cto;
  auto const ll = 2 * π * Radian * cto /
                  (primary_->angular_frequency() - nodal_precession_) / Day;

  Angle const ΔλE = -2 * π * Radian * cto / nto;
  Angle const δ = 2 * π * Radian / nto;

  int eto;
  for (int j = 1; j < cto; ++j) {
    if (mod(j * dto, cto) == 1 || mod(j * dto, -cto) == -1) {
      eto = j;
      break;
    }
  }

  for (int i = 0, k = 0; i < longitudes_of_ascending_nodes.size(); ++i) {
    Angle const Ω = longitudes_of_ascending_nodes[i];
    Instant const t = times_of_ascending_nodes[i];
    // TODO(egg): this assumes earthlike sign for the longitude.
    Angle const λ = Ω - (primary_->AngleAt(t) + π / 2 * Radian);
    if (i == 0) {
      k = geometry::Sign(ΔλE) * std::floor(λ / Abs(ΔλE));
    }
    terrestrial_longitudes_of_ascending_nodes.push_back(
        quantities::Mod(λ - (i + k) * ΔλE + π * Radian, 2 * π * Radian) -
        π * Radian);
  }

  std::vector<Angle> const λ =
      Unwind(terrestrial_longitudes_of_ascending_nodes);
  auto const λ0 = AverageOfCorrelated(λ);

  std::vector<Angle> λ_pe;
  std::vector<Angle> λ_ap;
  MeasurementResult<Angle> mean_λ_pe;
  MeasurementResult<Angle> mean_λ_ap;
  if (nto == 1 && cto == 1) {
    for (auto const& t : times_of_periapsides) {
      λ_pe.push_back(
          quantities::Mod((primary_centred_trajectory.EvaluatePosition(t) -
                           PrimaryCentred::origin)
                                  .coordinates()
                                  .ToSpherical()
                                  .longitude -
                              (primary_->AngleAt(t) + π / 2 * Radian) +
                              π * Radian,
                          2 * π * Radian) -
          π * Radian);
    }
    for (auto const& t : times_of_apoapsides) {
      λ_ap.push_back(
          quantities::Mod((primary_centred_trajectory.EvaluatePosition(t) -
                           PrimaryCentred::origin)
                                  .coordinates()
                                  .ToSpherical()
                                  .longitude -
                              (primary_->AngleAt(t) + π / 2 * Radian) +
                              π * Radian,
                          2 * π * Radian) -
          π * Radian);
    }
    mean_λ_pe = AverageOfCorrelated(Unwind(λ_pe));
    mean_λ_ap = AverageOfCorrelated(Unwind(λ_ap));
  }

  LOG(ERROR) << "--- General parameters ---";
  LOG(ERROR) << u8"T* = " << sidereal_period_ / Second << " s";
  LOG(ERROR) << u8"T* / T🜨 = " << sidereal_period_ / sidereal_rotation_period;
  LOG(ERROR) << u8"T🜨 / T* = " << sidereal_rotation_period / sidereal_period_;
  LOG(ERROR) << u8"T☊ = " << nodal_period_ / Second << " s";
  LOG(ERROR) << u8"T = " << anomalistic_period_ / Second << " s";
  LOG(ERROR) << u8"Ω′ = " << nodal_precession_ / (Degree / Day) << u8"°/d = "
             << nodal_precession_ / (Degree / JulianYear) << u8"°/a"
             << "\n = " << nodal_precession_ / (1 / tropical_year_) / Degree
             << u8"° / tropical year";
  LOG(ERROR) << u8"ω′ = " << apsidal_precession_ / (Degree / Day) << u8"°/d = "
             << apsidal_precession_ / (Degree / JulianYear) << u8"°/a";
  LOG(ERROR) << "i = " << inclination_ / Degree << u8"° (ψm)";
  LOG(ERROR) << "- Periods from osculating element regression -";
  LOG(ERROR) << "T_M = "
             << 2 * π * Radian /
                    LinearRegression(times, Unwind(mean_anomalies)).slope /
                    Second
             << " s";
  LOG(ERROR) << "T_l = "
             << 2 * π * Radian /
                    LinearRegression(times, Unwind(mean_longitudes)).slope /
                    Second
             << " s";
  LOG(ERROR)
      << "T_u = "
      << 2 * π * Radian /
             LinearRegression(times, Unwind(mean_arguments_of_latitude)).slope /
             Second
      << " s";
  LOG(ERROR) << "----";
  LOG(ERROR) << u8"  ± "
             << Variability(inclinations_at_extremal_latitudes,
                            inclination_.measured_value) /
                    Degree
             << u8"° (95%)";
  LOG(ERROR) << "e = " << eccentricity_;
  LOG(ERROR) << "e = " << 1 - ſ_r_pe_dt / ſ_a_dt
             << u8" (from integrated osculating a & r_pe)";
  LOG(ERROR) << "e = " << ſ_r_ap_dt / ſ_a_dt - 1
             << u8" (from integrated osculating a & r_ap)";
  LOG(ERROR) << "e = " << (ſ_r_ap_dt - ſ_r_pe_dt) / (ſ_r_ap_dt + ſ_r_pe_dt)
             << u8" (from integrated osculating r_pe & r_ap)";
  LOG(ERROR) << "  = " << (ſ_r_ap_dt - ſ_r_pe_dt) / (2 * ſ_a_dt)
             << u8" (from integrated osculating r_pe & r_ap, a = r_pe + r_ap)";
  LOG(ERROR) << u8"ω = "
             << quantities::ArcTan(ſ_e_sin_ω_dt, ſ_e_cos_ω_dt) / Degree
             << "° (from integrated eccentricity vector)";
  LOG(ERROR) << "e = "
             << quantities::Sqrt(
                    (Pow<2>(ſ_e_cos_ω_dt) + Pow<2>(ſ_e_sin_ω_dt))) /
                    (trajectory_.last().time() - trajectory_.Begin().time())
             << " (from integrated eccentricity vector)";
  LOG(ERROR) << "e = "
             << ſ_e_dt /
                    (trajectory_.last().time() - trajectory_.Begin().time())
             << " (from integrated osculating e)";
  LOG(ERROR) << "e = "
             << std::exp(ſ_log_e_dt / (trajectory_.last().time() -
                                       trajectory_.Begin().time()))
             << " (from integrated log osculating e)";
  LOG(ERROR) << "e = "
             << (trajectory_.last().time() - trajectory_.Begin().time()) /
                    ſ_e⁻¹_dt
             << " (from integrated inverse osculating e)";
  LOG(ERROR) << "a = "
             << ſ_a_dt /
                    (trajectory_.last().time() - trajectory_.Begin().time()) /
                    Kilo(Metre)
             << " km (integrated osculating)";
  LOG(ERROR) << "i = "
             << ſ_i_dt /
                    (trajectory_.last().time() - trajectory_.Begin().time()) /
                    Degree
             << u8"° (integrated osculating)";
  LOG(ERROR) << "--- Orbit with respect to the Earth ---";
  LOG(ERROR) << "- Phasing -";
  LOG(ERROR) << u8"κ = " << daily_recurrence_frequency;
  LOG(ERROR) << "N_To / C_To = " << nto << " / " << cto;
  LOG(ERROR) << u8"[ν0 ; DTo ; CTo] = [" << ν0 << " ; " << dto << " ; " << cto
             << "]";
  LOG(ERROR) << u8"𝕃 = NTo Td \t\t\t\t= " << nto * nodal_period_ / Day << " d";
  LOG(ERROR) << u8"𝕃 = 2π CTo / (Ω′T - Ω′) \t= " << ll << " d";

  LOG(ERROR) << u8"ΔλE = " << ΔλE / Degree << u8"°";
  LOG(ERROR) << u8"δ = " << δ / Degree << u8"°";
  LOG(ERROR) << u8"ETo* = " << eto;

  LOG(ERROR) << "- Ground track -";
  // TODOegg): extremal latitudes (& nominal inclination);
  LOG(ERROR) << u8"λ0 =" << λ0 / Degree << u8"°";
  LOG(ERROR) << u8"   ± " << Variability(λ, λ0.measured_value) / Degree
             << u8"° (95%)";
  LOG(ERROR) << u8"   ± "
             << Variability(λ, λ0.measured_value) *
                    ((OblateBody<Frame> const&)*primary_).reference_radius() /
                    Radian / Kilo(Metre)
             << u8" km (95%)";

  if (nto == 1 && cto == 1) {
    LOG(ERROR) << "- Geosynchronous orbit: central longitude -";
    LOG(ERROR) << u8"λ_pe =" << mean_λ_pe / Degree << u8"°";
    LOG(ERROR) << u8"   ± "
               << Variability(λ_pe, mean_λ_pe.measured_value) / Degree
               << u8"° (95%)";
    LOG(ERROR) << u8"λ_ap =" << mean_λ_ap / Degree << u8"°";
    LOG(ERROR) << u8"   ± "
               << Variability(λ_ap, mean_λ_ap.measured_value) / Degree
               << u8"° (95%)";
  }

  LOG(ERROR) << "- Orbit freezing -";
  LOG(ERROR) << u8"ω = " << ω / Degree << u8"°";
  LOG(ERROR) << u8"  ± "
             << Variability(arguments_of_periapsides, ω.measured_value) / Degree
             << u8"° (95%)";

  LOG(ERROR) << u8"r_p = " << periapsis_distance_ / Kilo(Metre) << " km";
  LOG(ERROR) << u8"    ± "
             << Variability(periapsis_distances,
                            periapsis_distance_.measured_value) /
                    Kilo(Metre)
             << u8" km (95%)";
  LOG(ERROR) << u8"r_a = " << apoapsis_distance_ / Kilo(Metre) << " km";
  LOG(ERROR) << u8"    ± "
             << Variability(apoapsis_distances,
                            apoapsis_distance_.measured_value) /
                    Kilo(Metre)
             << u8" km (95%)";
  LOG(ERROR) << u8"h_p = "
             << (periapsis_distance_ - primary_->mean_radius()) / Kilo(Metre)
             << " km";
  LOG(ERROR) << u8"h_a = "
             << (apoapsis_distance_ - primary_->mean_radius()) / Kilo(Metre)
             << " km";

  LOG(ERROR) << "--- Orbit with respect to the Sun, heliosynchronicity ---";

  // REMOVE BEFORE FLIGHT: There probably should be a sign here to turn the
  // tropical year into Ω'S.
  auto const ΔtS =
      2 * π * Radian /
      (nodal_precession_ - (2 * π * Radian / tropical_year_.measured_value));
  LOG(ERROR) << u8"ΔtS = " << ΔtS / Day << " d";

  LOG(ERROR) << u8"τNA = " << mean_τ / Degree << u8"° = "
             << 12 + (mean_τ * 24 / (2 * π * Radian)) << u8" h";
  LOG(ERROR) << u8"       ± "
             << Variability(τ, mean_τ.measured_value) * (24 / (2 * π * Radian))
             << " h (95 %)";
}

}  // namespace internal_orbit_analyser
}  // namespace astronomy
}  // namespace principia
