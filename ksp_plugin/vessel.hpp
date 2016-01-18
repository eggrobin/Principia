﻿
#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/part.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace integrators {
template <typename DifferentialEquation> class AdaptiveStepSizeIntegrator;
}  // namespace integrators
}  // namespace principia

namespace principia {

namespace ksp_plugin {
class Celestial;
template <typename InertialFrame, typename Frame> class Manœuvre;
}  // namespace ksp_plugin
namespace serialization {
class Vessel;
}  // namespace serialization

using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MasslessBody;
using quantities::GravitationalParameter;

namespace ksp_plugin {

// Represents a KSP |Vessel|.
class Vessel {
 public:
  using Manœuvres =
      std::vector<
          not_null<std::unique_ptr<Manœuvre<Barycentric, Navigation> const>>>;

  Vessel(Vessel const&) = delete;
  Vessel(Vessel&&) = delete;
  Vessel& operator=(Vessel const&) = delete;
  Vessel& operator=(Vessel&&) = delete;
  ~Vessel() = default;

  // Constructs a vessel whose parent is initially |*parent|.  No transfer of
  // ownership.
  explicit Vessel(not_null<Celestial const*> const parent);

  // Returns the body for this vessel.
  virtual not_null<MasslessBody const*> body() const;

  // True if, and only if, |history_| is not null.
  virtual bool is_synchronized() const;
  // True if, and only if, |prolongation_| is not null, i.e., if either
  // |CreateProlongation| or |CreateHistoryAndForkProlongation| was called at
  // some point.
  virtual bool is_initialized() const;

  virtual not_null<Celestial const*> parent() const;
  virtual void set_parent(not_null<Celestial const*> const parent);

  // Both accessors require |is_synchronized()|.
  virtual DiscreteTrajectory<Barycentric> const& history() const;
  virtual not_null<DiscreteTrajectory<Barycentric>*> mutable_history();

  // Both accessors require |is_initialized()|.
  virtual DiscreteTrajectory<Barycentric> const& prolongation() const;
  virtual not_null<DiscreteTrajectory<Barycentric>*> mutable_prolongation();

  // Requires |is_initialized()|.
  virtual not_null<FlightPlan*> flight_plan() const;
  virtual bool has_flight_plan() const;

  // Requires |has_prediction()|.
  virtual DiscreteTrajectory<Barycentric> const& prediction() const;
  virtual bool has_prediction() const;

  // Creates an |owned_prolongation_| for this vessel and appends a point with
  // the given |time| and |degrees_of_freedom|.  The vessel must not satisfy
  // |is_initialized()| nor |is_synchronized()|, |owned_prolongation_| must be
  // null.  The vessel |is_initialized()|, but does not satisfy
  // |is_synchronized()|, after the call.
  virtual void CreateProlongation(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Creates a |history_| for this vessel and appends a point with the
  // given |time| and |degrees_of_freedom|, then forks a |prolongation_| at
  // |time|.  Nulls |owned_prolongation_|.  The vessel must not satisfy
  // |is_synchronized()|.  |*owned_prolongation_| is destroyed *after*
  // |history_| has been constructed.
  // The vessel |is_synchronized()| and |is_initialized()| after the call.
  virtual void CreateHistoryAndForkProlongation(
      Instant const& time,
      DegreesOfFreedom<Barycentric> const& degrees_of_freedom);

  // Deletes the |prolongation_| and forks a new one at |time|.
  // The vessel must satisfy |is_synchronized()| and |is_initialized()|,
  // |owned_prolongation_| must be null.
  virtual void ResetProlongation(Instant const& time);

  // Creates a |flight_plan_| at the end of history using the given parameters.
  // Deletes any pre-existing predictions.
  // Does nothing unless |is_synchronized()|, pending the removal of
  // synchronization.
  // TODO(egg): struct containing (integrator, length tol, speed tol) so we
  // don't need that many parameters...
  virtual void CreateFlightPlan(
      Instant const& final_time,
      Mass const& initial_mass,
      not_null<Ephemeris<Barycentric>*> ephemeris,
      AdaptiveStepSizeIntegrator<
          Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance);

  // Deletes the |flight_plan_|.  Performs no action unless |has_flight_plan()|.
  virtual void DeleteFlightPlan();

  virtual void UpdatePrediction(
      not_null<Ephemeris<Barycentric>*> ephemeris,
      AdaptiveStepSizeIntegrator<
          Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
      Instant const& last_time,
      Length const& prediction_length_tolerance,
      Speed const& prediction_speed_tolerance);

  // Deletes the |prediction_|.  Performs no action unless |has_prediction()|.
  virtual void DeletePrediction();

  // The vessel must satisfy |is_initialized()|.
  virtual void WriteToMessage(
      not_null<serialization::Vessel*> const message) const;
  static not_null<std::unique_ptr<Vessel>> ReadFromMessage(
      serialization::Vessel const& message,
      not_null<Celestial const*> const parent);

 protected:
  // For mocking.
  Vessel();

 private:
  MasslessBody const body_;
  // The parent body for the 2-body approximation. Not owning.
  not_null<Celestial const*> parent_;
  // The past and present trajectory of the body. It ends at |HistoryTime()|
  // unless |*this| was created after |HistoryTime()|, in which case it ends
  // at |current_time_|.  It is advanced with a constant time step.
  std::unique_ptr<DiscreteTrajectory<Barycentric>> history_;
  // Most of the time, this is a child trajectory of |*history_|. It is forked
  // at |history_->last_time()| and continues until |current_time_|. It is
  // computed with a non-constant timestep, which breaks symplecticity.
  // If |history_| is null, this points to |owned_prolongation_| instead.
  // Not owning.
  DiscreteTrajectory<Barycentric>* prolongation_ = nullptr;
  // When the vessel is added, before it is synchonized with the other vessels
  // and celestials, there is no |history_|.  The prolongation is directly owned
  // during that time.  Null if, and only if, |history_| is not null.
  std::unique_ptr<DiscreteTrajectory<Barycentric>> owned_prolongation_;
  // Child trajectory of |history_|.
  DiscreteTrajectory<Barycentric>* prediction_ = nullptr;

  std::unique_ptr<FlightPlan> flight_plan_;
};

}  // namespace ksp_plugin
}  // namespace principia

#include "ksp_plugin/vessel_body.hpp"
