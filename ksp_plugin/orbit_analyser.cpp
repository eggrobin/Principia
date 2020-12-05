﻿
#include "ksp_plugin/orbit_analyser.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_orbit_analyser {

using base::dynamic_cast_not_null;
using base::MakeStoppableThread;
using geometry::Frame;
using geometry::NonRotating;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::DiscreteTrajectory;
using physics::KeplerOrbit;
using physics::MasslessBody;
using quantities::IsFinite;
using quantities::Infinity;

OrbitAnalyser::OrbitAnalyser(
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    Ephemeris<Barycentric>::FixedStepParameters analysed_trajectory_parameters)
    : ephemeris_(ephemeris),
      analysed_trajectory_parameters_(
          std::move(analysed_trajectory_parameters)) {}

void OrbitAnalyser::Restart() {
  analyser_ = MakeStoppableThread([this] { RepeatedlyAnalyseOrbit(); });
}

void OrbitAnalyser::RequestAnalysis(Parameters const& parameters) {
  if (!analyser_.joinable()) {
    analyser_ = MakeStoppableThread([this] { RepeatedlyAnalyseOrbit(); });
  }
  Ephemeris<Barycentric>::Guard guard(ephemeris_);
  if (ephemeris_->t_min() > parameters.first_time) {
    // Too much has been forgotten; we cannot perform this analysis.
    return;
  }
  last_parameters_ = parameters;
  absl::MutexLock l(&lock_);
  guarded_parameters_ = {std::move(guard), parameters};
}

std::optional<OrbitAnalyser::Parameters> const& OrbitAnalyser::last_parameters()
    const {
  return last_parameters_;
}

void OrbitAnalyser::RefreshAnalysis() {
  absl::MutexLock l(&lock_);
  if (next_analysis_.has_value()) {
    analysis_ = std::move(next_analysis_);
    next_analysis_.reset();
  }
}

OrbitAnalyser::Analysis* OrbitAnalyser::analysis() {
  return analysis_.has_value() ? &*analysis_ : nullptr;
}

double OrbitAnalyser::progress_of_next_analysis() const {
  return progress_of_next_analysis_;
}

Status OrbitAnalyser::RepeatedlyAnalyseOrbit() {
  for (;;) {
    // No point in going faster than 50 Hz.
    std::chrono::steady_clock::time_point const wakeup_time =
        std::chrono::steady_clock::now() + std::chrono::milliseconds(20);

    RETURN_IF_STOPPED;

    std::optional<GuardedParameters> guarded_parameters;
    {
      absl::MutexLock l(&lock_);
      if (!guarded_parameters_.has_value()) {
        // No parameters, let's wait for them to appear.
        continue;
      }
      std::swap(guarded_parameters, guarded_parameters_);
    }

    auto const& parameters = guarded_parameters->parameters;

    Analysis analysis{parameters.first_time};
    DiscreteTrajectory<Barycentric> trajectory;
    trajectory.Append(parameters.first_time,
                      parameters.first_degrees_of_freedom);

    RotatingBody<Barycentric> const* primary = nullptr;
    auto smallest_osculating_period = Infinity<Time>;
    for (auto const body : ephemeris_->bodies()) {
      RETURN_IF_STOPPED;
      auto const initial_osculating_elements =
          KeplerOrbit<Barycentric>{
              *body,
              MasslessBody{},
              parameters.first_degrees_of_freedom -
                  ephemeris_->trajectory(body)->EvaluateDegreesOfFreedom(
                      parameters.first_time),
              parameters.first_time}.elements_at_epoch();
      if (initial_osculating_elements.period.has_value() &&
          initial_osculating_elements.period < smallest_osculating_period) {
        smallest_osculating_period = *initial_osculating_elements.period;
        primary = dynamic_cast_not_null<RotatingBody<Barycentric> const*>(body);
      }
    }
    if (primary != nullptr) {
      std::vector<not_null<DiscreteTrajectory<Barycentric>*>> trajectories = {
          &trajectory};
      auto instance = ephemeris_->NewInstance(
          trajectories,
          Ephemeris<Barycentric>::NoIntrinsicAccelerations,
          analysed_trajectory_parameters_);
      Time const analysis_duration =
          std::min(parameters.extended_mission_duration.value_or(
                       parameters.mission_duration),
                   std::max(2 * smallest_osculating_period,
                            parameters.mission_duration));
      for (Instant t = parameters.first_time + analysis_duration / 0x1p10;
           trajectory.back().time < parameters.first_time + analysis_duration;
           t += analysis_duration / 0x1p10) {
        if (!ephemeris_->FlowWithFixedStep(t, *instance).ok()) {
          // TODO(egg): Report that the integration failed.
          break;
        }
        progress_of_next_analysis_ =
            (trajectory.back().time - parameters.first_time) /
            analysis_duration;
        RETURN_IF_STOPPED;
      }
      analysis.mission_duration_ =
          trajectory.back().time - parameters.first_time;

      // TODO(egg): |next_analysis_percentage_| only reflects the progress of
      // the integration, but the analysis itself can take a while; this results
      // in the progress bar being stuck at 100% while the elements and nodes
      // are being computed.

      using PrimaryCentred = Frame<enum class PrimaryCentredTag, NonRotating>;
      DiscreteTrajectory<PrimaryCentred> primary_centred_trajectory;
      BodyCentredNonRotatingDynamicFrame<Barycentric, PrimaryCentred>
          body_centred(ephemeris_, primary);
      for (auto const& [time, degrees_of_freedom] : trajectory) {
        RETURN_IF_STOPPED;
        primary_centred_trajectory.Append(
            time, body_centred.ToThisFrameAtTime(time)(degrees_of_freedom));
      }
      analysis.primary_ = primary;
      auto elements = OrbitalElements::ForTrajectory(
          primary_centred_trajectory, *primary, MasslessBody{});
      // We do not RETURN_IF_ERROR as ForTrajectory can return non-CANCELLED
      // statuses.
      RETURN_IF_STOPPED;
      if (elements.ok()) {
        analysis.elements_ = std::move(elements).ValueOrDie();
        // TODO(egg): max_abs_Cᴛₒ should probably depend on the number of
        // revolutions.
        auto recurrence = OrbitRecurrence::ClosestRecurrence(
            analysis.elements_->nodal_period(),
            analysis.elements_->nodal_precession(),
            *primary,
            /*max_abs_Cᴛₒ=*/100);
        if (recurrence.ok()) {
          analysis.closest_recurrence_ = std::move(recurrence).ValueOrDie();
        }
        auto ground_track =
            OrbitGroundTrack::ForTrajectory(primary_centred_trajectory,
                                            *primary,
                                            /*mean_sun=*/std::nullopt);
        RETURN_IF_ERROR(ground_track);
        analysis.ground_track_ = std::move(ground_track).ValueOrDie();
        analysis.ResetRecurrence();
      }
    }

    {
      absl::MutexLock l(&lock_);
      next_analysis_ = std::move(analysis);
    }

    std::this_thread::sleep_until(wakeup_time);
  }
}

Instant const& OrbitAnalyser::Analysis::first_time() const {
  return first_time_;
}

Time const& OrbitAnalyser::Analysis::mission_duration() const {
  return mission_duration_;
}

RotatingBody<Barycentric> const* OrbitAnalyser::Analysis::primary() const {
  return primary_;
}

std::optional<OrbitalElements> const& OrbitAnalyser::Analysis::elements()
    const {
  return elements_;
}

std::optional<OrbitRecurrence> const& OrbitAnalyser::Analysis::recurrence()
    const {
  return recurrence_;
}

std::optional<OrbitGroundTrack> const& OrbitAnalyser::Analysis::ground_track()
    const {
  return ground_track_;
}

std::optional<OrbitGroundTrack::EquatorCrossingLongitudes> const&
OrbitAnalyser::Analysis::equatorial_crossings() const {
  return equatorial_crossings_;
}

void OrbitAnalyser::Analysis::SetRecurrence(
    OrbitRecurrence const& recurrence) {
  if (recurrence_ != recurrence) {
    recurrence_ = recurrence;
    if (ground_track_.has_value()) {
      equatorial_crossings_ = ground_track_->equator_crossing_longitudes(
          recurrence, /*first_ascending_pass_index=*/1);
    }
  }
}

void OrbitAnalyser::Analysis::ResetRecurrence() {
  if (closest_recurrence_.has_value()) {
    SetRecurrence(*closest_recurrence_);
  } else {
    recurrence_.reset();
    equatorial_crossings_.reset();
  }
}

OrbitAnalyser::Analysis::Analysis(Instant const& first_time)
    : first_time_(first_time) {}

}  // namespace internal_orbit_analyser
}  // namespace ksp_plugin
}  // namespace principia
