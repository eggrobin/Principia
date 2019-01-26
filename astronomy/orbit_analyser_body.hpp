#pragma once

#include "astronomy/orbit_analyser.hpp"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "physics/apsides.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/kepler_orbit.hpp"
#include "testing_utilities/statistics.hpp"

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
using quantities::Length;
using quantities::Speed;
using quantities::Time;
using quantities::si::Radian;
using testing_utilities::Mean;
using testing_utilities::Slope;

// TODO(egg): move this to a new file quantities/uncertainty.hpp.
using quantities::NaN;
using quantities::Difference;
using quantities::Square;
using quantities::Sqrt;

template<typename T>
struct MeasurementResult {
  T measured_value;
  Difference<T> standard_measurement_uncertainty;
};

template<typename T>
T Average(std::vector<T> const& samples) {
  // TODO(egg): BarycentreCalculator should have an unweighted version.
  T total{};
  for (T const& sample : samples) {
    total += sample - T{};
  }
  return total / samples.size();
}

// See Guide to the expression of uncertainty in measurement, 4.2.
template<typename T>
MeasurementResult<T> AverageOfIndependent(std::vector<T> const& measured_values) {
  int const n = measured_values.size();
  T average = Average(measured_values);
  Square<Difference<T>> experimental_variance{};
  for (T const& measured_value : measured_values) {
    experimental_variance += Pow<2>(measured_value - average);
  }
  experimental_variance /= n - 1;
  Difference<T> experimental_standard_deviation_of_the_mean =
      Sqrt(experimental_variance / n);
  return {average, experimental_standard_deviation_of_the_mean};
}

template<typename T>
MeasurementResult<T> AverageOfCorrelated(std::vector<T> const& time_series) {
  int const n = time_series.size();
  T average = Average(time_series);
  // See:
  // — Grant Foster (1996), Time Series Analysis by Projection. I. Statistical
  //   Properties of Fourier Analysis, equation 5.4;
  // — Chris Koen and Fred Lombard (1993), The analysis of indexed astronomical
  //   time series — I. Basic methods, equations 15, 19, and 2 (we use L = 1).
  // Note that there is a typo in equation 2 of Koen & Lombard, the t is missing
  // from the exponent.  A correct expression for the periodogram may be found
  // in equation 4 of Chris Koen and Fred Lombard (2004), The analysis of
  // indexed astronomical time series — IX. A period change test.
  Difference<T> re_second_fourier_coefficient{};
  Difference<T> im_second_fourier_coefficient{};
  // We use the convention from Koen and Lombard (2004), ω = 2π/n, rather than
  // ω = 1/n as in Koen and Lombard (1993).
  double const ω = 2 * π / n;
  for (int t = 0; t < n; ++t) {
    Difference<T> const Δy = time_series[t] - average;
    re_second_fourier_coefficient += Δy * std::cos(ω * t);
    im_second_fourier_coefficient += Δy * std::sin(ω * t);
  }
  Square<Difference<T>> const spectral_density_at_0_frequency =
      (Pow<2>(re_second_fourier_coefficient) +
       Pow<2>(im_second_fourier_coefficient)) / n;
  return {average, Sqrt(spectral_density_at_0_frequency / n)};
}

std::ostream& operator<<(std::ostream& out,
                         MeasurementResult<double> measurement) {
  int const value_decimal_exponent =
      std::floor(std::log10(measurement.measured_value));
  int const uncertainty_decimal_exponent =
      std::floor(std::log10(measurement.standard_measurement_uncertainty));
  int const significant_digits =
      value_decimal_exponent - uncertainty_decimal_exponent;
  std::string value_digits = std::to_string(static_cast<int>(std::nearbyint(
      std::abs(measurement.measured_value) *
      std::pow(10, (significant_digits + 2) - value_decimal_exponent - 1))));
  std::string uncertainty_digits = std::to_string(static_cast<int>(
      std::nearbyint(std::abs(measurement.standard_measurement_uncertainty) *
                     std::pow(10, 2 - uncertainty_decimal_exponent - 1))));
  return out << (std::signbit(measurement.measured_value) ? "-" : "+")
             << value_digits[0] << "." << value_digits.substr(1) << "("
             << uncertainty_digits << u8") × 10^"
             << (value_decimal_exponent >= 0 ? "+" : "")
             << value_decimal_exponent;
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

  // Try to make things scale-free by deriving tolerances from the orbit.
  // TODO(egg): perhaps we want to use a fixed-step integrator instead?
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
      t0 + Abs(2 * π * Radian / nodal_precession_),
      parameters,
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/false);
  RecomputeProperties();

  ephemeris_->FlowWithAdaptiveStep(
      &trajectory_,
      Ephemeris<Frame>::NoIntrinsicAcceleration,
      t0 + Abs(2 * π * Radian / apsidal_precession_),
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

  // TODO(egg): the statistics functions should support affine spaces:
  // |times_of_ascending_nodes| should hold instants.
  std::vector<Time> times_of_periapsides;
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
      times_between_periapsides.push_back(periapsis.time() - J2000 -
                                          times_of_periapsides.back());
    }
    times_of_periapsides.push_back(periapsis.time() - J2000);
    arguments_of_periapsides.push_back(ω);
  }
  apsidal_precession_ = Slope(times_of_periapsides, arguments_of_periapsides);
  anomalistic_period_ = Mean(times_between_periapsides);
  LOG(ERROR) << u8"ω′ = " << apsidal_precession_;
  LOG(ERROR) << u8"T = " << anomalistic_period_;
  LOG(ERROR) << "n = " << times_between_periapsides.size();
  auto T = AverageOfCorrelated(times_between_periapsides);
  LOG(ERROR) << u8"T = "
             << MeasurementResult<double>{T.measured_value / Second,
                                          T.standard_measurement_uncertainty /
                                              Second}
             << " s";
  base::OFStream f(SOLUTION_DIR / ("anomalistic_period_" + std::to_string(times_between_periapsides.size())));
  f << mathematica::Assign("tPe" + std::to_string(times_between_periapsides.size()), times_between_periapsides); 


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

  // TODO(egg): the statistics functions should support affine spaces:
  // |times_of_ascending_nodes| should hold instants.
  std::vector<Time> times_of_ascending_nodes;
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
      times_between_ascending_nodes.push_back(node.time() - J2000 -
                                              times_of_ascending_nodes.back());
    }
    times_of_ascending_nodes.push_back(node.time() - J2000);
    longitudes_of_ascending_nodes.push_back(Ω);
  }
  nodal_precession_ =
      Slope(times_of_ascending_nodes, longitudes_of_ascending_nodes);
  nodal_period_ = Mean(times_between_ascending_nodes);
  LOG(ERROR) << u8"Ω′ = " << nodal_precession_;
  LOG(ERROR) << u8"T☊ = " << nodal_period_;

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
  sidereal_period_ = Mean(times_between_xz_ascensions);
  LOG(ERROR) << u8"T* = " << sidereal_period_;
  Time const sidereal_rotation_period =
      2 * π * Radian / primary_->angular_frequency();
  LOG(ERROR) << u8"T* / T🜨 = " << sidereal_period_ / sidereal_rotation_period;
  LOG(ERROR) << u8"T🜨 / T* = " << sidereal_rotation_period / sidereal_period_;
}

}  // namespace internal_orbit_analyser
}  // namespace astronomy
}  // namespace principia
