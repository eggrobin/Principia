#pragma once

#include "astronomy/orbit_analyser.hpp"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "physics/apsides.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/kepler_orbit.hpp"

namespace principia {
namespace astronomy {
namespace internal_orbit_analyser {

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
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Radian;
using quantities::si::Second;

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
      *initial_osculating_elements.semimajor_axis * 1e-6;
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
      t0 + 10 * initial_osculating_period,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 100 * initial_osculating_period,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + 1000 * initial_osculating_period,
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + Abs(2 * π * Radian / nodal_precession_.measured_value),
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + Abs(2 * π * Radian / apsidal_precession_.measured_value),
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();
}

template<typename Frame>
void OrbitAnalyser<Frame>::RecomputeProperties() {

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
  for (auto periapsis = periapsides.Begin();
       periapsis != periapsides.End();
       ++periapsis) {
    // TODO(egg): We could probably do something more efficient, because we know
    // that we are at the periapsis, and we only need the argument of periapsis.
    Angle ω = *KeplerOrbit<Frame>(
                  *primary_,
                  MasslessBody{},
                  ephemeris_->trajectory(primary_)->EvaluateDegreesOfFreedom(
                      periapsis.time()) -
                      periapsis.degrees_of_freedom(),
                  periapsis.time()).elements_at_epoch().argument_of_periapsis;
    if (!arguments_of_periapsides.empty()) {
      ω += std::nearbyint((arguments_of_periapsides.back() - ω) /
                          (2 * π * Radian)) *
           2 * π * Radian;
      times_between_periapsides.push_back(periapsis.time() -
                                          times_of_periapsides.back());
    }
    times_of_periapsides.push_back(periapsis.time());
    arguments_of_periapsides.push_back(ω);
  }
  apsidal_precession_ =
      LinearRegression(times_of_periapsides, arguments_of_periapsides).slope;
  anomalistic_period_ = AverageOfCorrelated(times_between_periapsides);
  LOG(ERROR) << u8"ω′ = " << apsidal_precession_ / (Degree / Day) << u8"° / d";
  LOG(ERROR) << u8"T = " << anomalistic_period_ / Second << " s";
  LOG(ERROR) << "n = " << times_between_periapsides.size();

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

  DiscreteTrajectory<PrimaryCentred> ascending_nodes;
  DiscreteTrajectory<PrimaryCentred> descending_nodes;
  ComputeNodes(primary_centred_trajectory.Begin(),
               primary_centred_trajectory.End(),
               Vector<double, PrimaryCentred>({0, 0, 1}),
               ascending_nodes,
               descending_nodes);

  std::vector<Instant> times_of_ascending_nodes;
  std::vector<Time> times_between_ascending_nodes;
  std::vector<Angle> longitudes_of_ascending_nodes;
  for (auto node = ascending_nodes.Begin(); node != ascending_nodes.End(); ++node) {
    // We do not construct |KeplerianElements|: we only need the longitude of
    // the ascending node, and we are at the ascending node so the computation
    // is trivial.
    // In order to compute a linear fit, we have to add or subtract
    // cycles as appropriate as the angle varies.
    Angle Ω = OrientedAngleBetween(
        x, node.degrees_of_freedom().position() - PrimaryCentred::origin, z);
    if (!longitudes_of_ascending_nodes.empty()) {
      Ω += std::nearbyint((longitudes_of_ascending_nodes.back() - Ω) /
                          (2 * π * Radian)) *
           2 * π * Radian;
      times_between_ascending_nodes.push_back(node.time() -
                                              times_of_ascending_nodes.back());
    }
    times_of_ascending_nodes.push_back(node.time());
    longitudes_of_ascending_nodes.push_back(Ω);
  }
  nodal_precession_ =
      LinearRegression(times_of_ascending_nodes,
                       longitudes_of_ascending_nodes).slope;
  nodal_period_ = AverageOfCorrelated(times_between_ascending_nodes);
  LOG(ERROR) << u8"Ω′ = " << nodal_precession_ / (Degree / Day) << u8"° / d";;
  LOG(ERROR) << u8"T☊ = " << nodal_period_ / Second << " s";

  // By computing the "nodes" with respect to the xz plane, i.e., the crossings
  // of the xz plane, we compute the points at which the projection of the orbit
  // onto the reference plane xy passes the fixed reference direction x.
  // While these points themselves are largely meaningless, the time between
  // them is the sidereal period.
  DiscreteTrajectory<PrimaryCentred> xz_ascensions;
  DiscreteTrajectory<PrimaryCentred> xz_descents;
  ComputeNodes(primary_centred_trajectory.Begin(),
               primary_centred_trajectory.End(),
               z * x,
               xz_ascensions,
               xz_descents);
  std::vector<Time> times_between_xz_ascensions;
  auto ascension = xz_ascensions.Begin();
  Instant previous_time = ascension.time();
  for (; ascension != xz_ascensions.End();
       previous_time = ascension.time(), ++ascension) {
    times_between_xz_ascensions.push_back(ascension.time() - previous_time);
  }
  sidereal_period_ = AverageOfCorrelated(times_between_xz_ascensions);
  LOG(ERROR) << u8"T* = " << sidereal_period_ / Second << " s";
  Time const sidereal_rotation_period =
      2 * π * Radian / primary_->angular_frequency();
  LOG(ERROR) << u8"T* / T🜨 = " << sidereal_period_ / sidereal_rotation_period;
  LOG(ERROR) << u8"T🜨 / T* = " << sidereal_rotation_period / sidereal_period_;
}

}  // namespace internal_orbit_analyser
}  // namespace astronomy
}  // namespace principia
