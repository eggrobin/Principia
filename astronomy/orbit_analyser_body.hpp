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
using geometry::OrientedAngleBetween;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using geometry::Bivector;
using geometry::Velocity;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::ComputeApsides;
using physics::ComputeNodes;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MasslessBody;
using physics::RigidMotion;
using physics::RigidTransformation;
using quantities::Abs;
using quantities::Angle;
using quantities::AverageOfCorrelated;
using quantities::Length;
using quantities::LinearRegression;
using quantities::Speed;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Radian;
using quantities::si::Second;

// |angles| should be sampled from a slowly-varying continuous function
// f: ℝ →  𝑆¹ = ℝ / 2π ℝ (specifically, consecutive angles should  differ by
// less than π).  Returns the corresponding sampling of the continuous g: ℝ → ℝ
// such that f = g mod 2π and f(0) = g(0).
std::vector<Angle> Unwind(std::vector<Angle> const& angles) {
  if (angles.empty()) {
    return angles;
  }
  std::vector<Angle> unwound_angles;
  unwound_angles.reserve(angles.size());
  unwound_angles.push_back(angles.front());
  for (int i = 1; i < angles.size(); ++i) {
    unwound_angles.push_back(
        angles[i] +
        std::nearbyint((unwound_angles.back() - angles[i]) / (2 * π * Radian)) *
            2 * π * Radian);
  }
  return unwound_angles;
}

template<typename Frame>
OrbitAnalyser<Frame>::OrbitAnalyser(
    not_null<Ephemeris<Frame>*> const ephemeris,
    not_null<RotatingBody<Frame> const*> const primary,
    Instant const initial_time,
    DegreesOfFreedom<Frame> const initial_state)
    : ephemeris_(ephemeris), primary_(primary) {
  trajectory_.Append(initial_time, initial_state);
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

LOG(ERROR) << 10;
  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 10 * initial_osculating_period,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

LOG(ERROR) << 100;
  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 100 * nodal_period_.measured_value,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

LOG(ERROR) << 1000;
  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 1000 * nodal_period_.measured_value,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

LOG(ERROR) << "1 nodal precession cycle";
  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + Abs(2 * π * Radian / nodal_precession_.measured_value),
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();


LOG(ERROR) << "10'000";
  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 10'000 * nodal_period_.measured_value,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

LOG(ERROR) << "20'000";
  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 20'000 * nodal_period_.measured_value,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

LOG(ERROR) << "30'000";
  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 30'000 * nodal_period_.measured_value,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();
}

template<typename Frame>
void OrbitAnalyser<Frame>::RecomputeProperties() {
  LOG(ERROR) << (trajectory_.last().time() - trajectory_.Begin().time()) / Day
             << " d";
  LOG(ERROR) << (trajectory_.last().time() - trajectory_.Begin().time()) /
                    JulianYear
             << " a";

  DiscreteTrajectory<Frame> periapsides;
  DiscreteTrajectory<Frame> apoapsides;
  ComputeApsides(*ephemeris_->trajectory(primary_),
                 trajectory_.Begin(),
                 trajectory_.End(),
                 periapsides,
                 apoapsides);

  std::vector<Instant> times_of_periapsides;
  std::vector<Time> times_between_periapsides;
  std::vector<Angle> arguments_of_periapsides;
  std::vector<Length> periapsis_distances;
  std::vector<Length> apoapsis_distances;

  for (auto periapsis = periapsides.Begin();
       periapsis != periapsides.End();
       ++periapsis) {
    // TODO(egg): We could probably do something a lot more efficient, because
    // we know that we are at the periapsis, and we only need the argument of
    // periapsis.
    Angle ω = *KeplerOrbit<Frame>(
                  *primary_,
                  MasslessBody{},
                  ephemeris_->trajectory(primary_)->EvaluateDegreesOfFreedom(
                      periapsis.time()) -
                      periapsis.degrees_of_freedom(),
                  periapsis.time()).elements_at_epoch().argument_of_periapsis;
    if (!times_of_periapsides.empty()) {
      times_between_periapsides.push_back(periapsis.time() -
                                          times_of_periapsides.back());
    }
    times_of_periapsides.push_back(periapsis.time());
    arguments_of_periapsides.push_back(ω);
    periapsis_distances.push_back(
        (ephemeris_->trajectory(primary_)->EvaluatePosition(periapsis.time()) -
         periapsis.degrees_of_freedom().position()).Norm());
  }
  for (auto apoapsis = apoapsides.Begin();
       apoapsis != apoapsides.End();
       ++apoapsis) {
    apoapsis_distances.push_back(
        (ephemeris_->trajectory(primary_)->EvaluatePosition(apoapsis.time()) -
         apoapsis.degrees_of_freedom().position()).Norm());
  }
  apsidal_precession_ =
      LinearRegression(times_of_periapsides,
                       Unwind(arguments_of_periapsides)).slope;
  anomalistic_period_ = AverageOfCorrelated(times_between_periapsides);
  periapsis_distance_ = AverageOfCorrelated(periapsis_distances);
  apoapsis_distance_ = AverageOfCorrelated(apoapsis_distances);
  LOG(ERROR) << u8"ω′ = " << apsidal_precession_ / (Degree / JulianYear)
             << u8"°/a";
  LOG(ERROR) << u8"T = " << anomalistic_period_ / Second << " s";
  LOG(ERROR) << u8"r_p = " << periapsis_distance_ / Kilo(Metre) << " km";
  LOG(ERROR) << u8"r_a = " << apoapsis_distance_ / Kilo(Metre) << " km";
  LOG(ERROR) << u8"h_p = " << (periapsis_distance_ - primary_->mean_radius()) / Kilo(Metre) << " km";
  LOG(ERROR) << u8"h_a = " << (apoapsis_distance_ - primary_->mean_radius()) / Kilo(Metre) << " km";

  enum class PrimaryTag { normal, sideways };
  // The origin of the reference frame is the centre of mass of |*primary_|.
  // The axes are those of |Frame|.
  // The reference plane for orbital analysis is the xy plane.
  using PrimaryCentred = geometry::Frame<PrimaryTag,
                                         PrimaryTag::normal,
                                         /*frame_is_inertial=*/false>;
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
          Rotation<PrimaryCentred, PrimaryCentredSideways>(
              x, z, x * z).Forget()),
      /*angular_velocity_of_to_frame=*/AngularVelocity<PrimaryCentred>{},
      /*velocity_of_to_frame_origin=*/
      Velocity<PrimaryCentred>{});

  BodyCentredNonRotatingDynamicFrame<Frame, PrimaryCentred> primary_centred(
      ephemeris_, primary_);

  DiscreteTrajectory<PrimaryCentred> primary_centred_trajectory;
  DiscreteTrajectory<PrimaryCentredSideways> sideways_primary_centred_trajectory;
  for (auto it = trajectory_.Begin(); it != trajectory_.End(); ++it) {
    primary_centred_trajectory.Append(
        it.time(),
        primary_centred.ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
    sideways_primary_centred_trajectory.Append(
        it.time(), tip(primary_centred_trajectory.last().degrees_of_freedom()));
  }

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
  base::OFStream f(SOLUTION_DIR / "sidereal_period_fine");
  f << mathematica::Assign("t", times_between_xz_ascensions);
  sidereal_period_ = AverageOfCorrelated(times_between_xz_ascensions);
  LOG(ERROR) << u8"T* = " << sidereal_period_ / Second << " s";
  Time const sidereal_rotation_period =
      2 * π * Radian / primary_->angular_frequency();
  LOG(ERROR) << u8"T* / T🜨 = " << sidereal_period_ / sidereal_rotation_period;
  LOG(ERROR) << u8"T🜨 / T* = " << sidereal_rotation_period / sidereal_period_;

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
  std::vector<Angle> terrestrial_longitudes_of_ascending_nodes;
  for (auto node = ascending_nodes.Begin(); node != ascending_nodes.End(); ++node) {
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
  }
  nodal_precession_ =
      LinearRegression(times_of_ascending_nodes,
                       Unwind(longitudes_of_ascending_nodes)).slope;
  nodal_period_ = AverageOfCorrelated(times_between_ascending_nodes);

  LOG(ERROR) << u8"Ω′ = " << nodal_precession_ / (Degree / JulianYear)
             << u8"°/a";
  LOG(ERROR) << u8"T☊ = " << nodal_period_ / Second << " s";

  // (7.41).
  MeasurementResult<double> const daily_recurrence_frequency =
      (2 * π * Radian / nodal_period_) /
      (primary_->angular_frequency() - nodal_precession_);
  LOG(ERROR) << u8"κ = " << daily_recurrence_frequency;

  // 11.7.2.
  double smallest_fraction = std::numeric_limits<double>::infinity();
  int cycle_days;
  int nto;
  int cto;
  for (int j = 1; j < 50; ++j) {
    MeasurementResult<double> κ_j = daily_recurrence_frequency * j;
    if (κ_j.standard_uncertainty > 0.5) {
      LOG(ERROR) << "within uncertainty at J = " << j;
    }
    double const abs_κ_j = std::abs(κ_j.measured_value);
    double const fraction = std::abs(abs_κ_j - std::nearbyint(abs_κ_j));
    if (fraction < smallest_fraction) {
      cycle_days = j;
      smallest_fraction = fraction;
      nto = std::nearbyint(abs_κ_j);
      cto = j;
      LOG(ERROR) << u8"frac |κJ| = " << fraction;
      LOG(ERROR) << "for J = " << j;
    }
  }
  LOG(ERROR) << "N_To = " << nto;
  LOG(ERROR) << "C_To = " << cto;

  int const ν0 = std::nearbyint(daily_recurrence_frequency.measured_value);
  int const dto = nto - ν0 * cto;
  LOG(ERROR) << u8"[ν0 ; DTo ; CTo] = ["
             << ν0 << " ; " << dto << " ; " << cto << "]";
  LOG(ERROR) << u8"𝕃 = NTo Td \t\t\t\t= " << nto * nodal_period_ / Day << " d";
  auto const ll = 2 * π * Radian * cto /
                    (primary_->angular_frequency() - nodal_precession_) / Day;
  LOG(ERROR) << u8"𝕃 = 2π CTo / (Ω′T - Ω′) \t= "
             << ll 
             << " d";

  Angle const ΔλE = -2 * π * Radian * cto / nto;
  Angle const δ = 2 * π * Radian / nto;

  int eto;
  for (int j = 1; j < cto; ++j) {
    if (mod(j * dto, cto) == 1 || mod(j * dto, -cto) == -1) {
      eto = j;
      break;
    }
  }

  LOG(ERROR) << u8"ΔλE = " << ΔλE / Degree << u8"°";
  LOG(ERROR) << u8"δ = " << δ / Degree << u8"°";
  LOG(ERROR) << u8"ETo* = " << eto;

  for (int i = 0; i < longitudes_of_ascending_nodes.size(); ++i) {
    Angle const Ω = longitudes_of_ascending_nodes[i];
    Instant const t = times_of_ascending_nodes[i];
    // TODO(egg): this assumes earthlike sign for the longitude.
    Angle const λ = Ω - (primary_->AngleAt(t) + π / 2 * Radian) + π * Radian;
    terrestrial_longitudes_of_ascending_nodes.push_back(
        quantities::Mod(
            λ - i * ΔλE,
            2 * π * Radian) -
        π * Radian);
  }

  std::vector<Angle> const λ = Unwind(terrestrial_longitudes_of_ascending_nodes);
  auto const λ0 = AverageOfCorrelated(λ);
  auto const λmax = *std::max_element(λ.begin(), λ.end());
  auto const λmin = *std::min_element(λ.begin(), λ.end());
  LOG(ERROR) << u8"λ0 =" << λ0 / Degree << u8"°";
  LOG(ERROR) << u8"λ- =" << λmin / Degree << u8"°";
  LOG(ERROR) << u8"λ+ =" << λmax / Degree << u8"°";
  LOG(ERROR) << u8"Δλ =" << (λmax - λmin) / Degree << u8"°";
  LOG(ERROR) << mathematica::ToMathematica(λ);
  base::OFStream tf(SOLUTION_DIR / "longitudes");
  tf << mathematica::Assign("longitudes", terrestrial_longitudes_of_ascending_nodes);
  tf << mathematica::Assign("t", times_between_xz_ascensions);
}

}  // namespace internal_orbit_analyser
}  // namespace astronomy
}  // namespace principia
