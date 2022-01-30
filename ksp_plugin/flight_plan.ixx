module;

#include <vector>

#include "absl/status/status.h"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/discrete_trajectory_segment_iterator.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

export module principia.ksp_plugin.flight_plan;

import principia.ksp_plugin.frames;
import principia.ksp_plugin.man�uvre;
import principia.ksp_plugin.orbit_analyser;

namespace principia {

using base::not_null;
using geometry::Instant;
using integrators::AdaptiveStepSizeIntegrator;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::DiscreteTrajectorySegmentIterator;
using physics::Ephemeris;
using quantities::Length;
using quantities::Mass;
using quantities::Speed;

export namespace ksp_plugin {

// A chain of trajectories obtained by executing the corresponding
// |NavigationMan�uvre|s.
class FlightPlan {
 public:
  // Creates a |FlightPlan| with no burns starting at |initial_time| with
  // |initial_degrees_of_freedom| and with the given |initial_mass|.  The
  // trajectories are computed using the given parameters by the given
  // |ephemeris|.  The flight plan contains a single coast which, if possible
  // ends at |desired_final_time|.
  FlightPlan(Mass const& initial_mass,
             Instant const& initial_time,
             DegreesOfFreedom<Barycentric> initial_degrees_of_freedom,
             Instant const& desired_final_time,
             not_null<Ephemeris<Barycentric>*> ephemeris,
             Ephemeris<Barycentric>::AdaptiveStepParameters
                 adaptive_step_parameters,
             Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
                 generalized_adaptive_step_parameters);
  virtual ~FlightPlan() = default;

  // Construction parameters.
  virtual Instant initial_time() const;
  virtual Instant desired_final_time() const;

  // End time of the last coast.  If this is less than |desired_final_time()|,
  // there is at least an anomalous man�uvre.
  virtual Instant actual_final_time() const;

  // The number of man�uvres in the flight plan.
  virtual int number_of_man�uvres() const;

  // The number of man�uvres at the end of flight plan that are anomalous, i.e.,
  // follow an anomalous segment.  These are man�uvres whose Frenet trihedron
  // cannot be drawn.  The functions that change man�uvres may change the number
  // of anomalous man�uvres.
  virtual int number_of_anomalous_man�uvres() const;

  // Returns the specified man�uvre.  |index| must be in
  // [0, number_of_man�uvres()[.
  virtual NavigationMan�uvre const& GetMan�uvre(int index) const;

  // Inserts a man�uvre at the given |index| using the specified |burn|. |index|
  // must be in [0, number_of_man�uvres()].  Returns an error and has no effect
  // if the given |burn| cannot fit between the preceding burn and the following
  // one or the end of the flight plan.  Otherwise, updates the flight plan and
  // returns the integration status.
  virtual absl::Status Insert(NavigationMan�uvre::Burn const& burn, int index);

  // Removes the man�uvre with the given |index|, which must be in
  // [0, number_of_man�uvres()[.
  virtual absl::Status Remove(int index);

  // Replaces a man�uvre with one using the specified |burn|.  |index| must be
  // in [0, number_of_man�uvres()[.  Returns an error and has no effect if the
  // given |burn| cannot fit between the preceding and following burns.
  // Otherwise, updates the flight plan and returns the integration status.
  virtual absl::Status Replace(NavigationMan�uvre::Burn const& burn, int index);

  // Updates the desired final time of the flight plan.  Returns an error and
  // has no effect if |desired_final_time| is before the beginning of the last
  // coast.
  virtual absl::Status SetDesiredFinalTime(Instant const& desired_final_time);

  // Sets the parameters used to compute the trajectories and recomputes them
  // all.  Returns the integration status.
  virtual absl::Status SetAdaptiveStepParameters(
      Ephemeris<Barycentric>::AdaptiveStepParameters const&
          adaptive_step_parameters,
      Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
          generalized_adaptive_step_parameters);

  virtual Ephemeris<Barycentric>::AdaptiveStepParameters const&
  adaptive_step_parameters() const;
  virtual Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
  generalized_adaptive_step_parameters() const;

  // Returns the number of trajectories in this object.
  virtual int number_of_segments() const;

  // |index| must be in [0, number_of_segments()[.
  virtual DiscreteTrajectorySegmentIterator<Barycentric>
  GetSegment(int index) const;
  virtual DiscreteTrajectory<Barycentric> const& GetAllSegments() const;

  // |coast_index| must be in [0, number_of_man�uvres()].
  virtual OrbitAnalyser::Analysis* analysis(int coast_index);
  double progress_of_analysis(int coast_index) const;

  void WriteToMessage(not_null<serialization::FlightPlan*> message) const;

  // This may return a null pointer if the flight plan contained in the
  // |message| is anomalous.
  static std::unique_ptr<FlightPlan> ReadFromMessage(
      serialization::FlightPlan const& message,
      not_null<Ephemeris<Barycentric>*> ephemeris);

  static constexpr std::int64_t max_ephemeris_steps_per_frame = 1000;

  static constexpr absl::StatusCode bad_desired_final_time =
      absl::StatusCode::kOutOfRange;
  static constexpr absl::StatusCode does_not_fit =
      absl::StatusCode::kOutOfRange;
  static constexpr absl::StatusCode singular =
      absl::StatusCode::kInvalidArgument;

 protected:
  // For mocking.
  FlightPlan();

 private:
  // Clears and recomputes all trajectories in |segments_|.
  absl::Status RecomputeAllSegments();

  // Flows the given |segment| for the duration of |man�uvre| using its
  // intrinsic acceleration.
  absl::Status BurnSegment(
      NavigationMan�uvre const& man�uvre,
      DiscreteTrajectorySegmentIterator<Barycentric> segment);

  // Flows the given |segment| until |desired_final_time| with no intrinsic
  // acceleration.
  absl::Status CoastSegment(
      Instant const& desired_final_time,
      DiscreteTrajectorySegmentIterator<Barycentric> segment);

  // Computes new trajectories and appends them to |segments_|.  This updates
  // the last coast of |segments_| and then appends one coast and one burn for
  // each man�uvre in |man�uvres|.  If one of the integration returns an error,
  // returns that error.  In this case the trajectories that follow the one in
  // error are of length 0 and are anomalous.
  // TODO(phl): The argument should really be an std::span, but then Apple has
  // invented the Macintosh.
  absl::Status ComputeSegments(std::vector<NavigationMan�uvre>::iterator begin,
                               std::vector<NavigationMan�uvre>::iterator end);

  // Adds a trajectory to |segments_|, forked at the end of the last one.  If
  // there are already anomalous trajectories, the newly created trajectory is
  // anomalous too.
  void AddLastSegment();

  // Forgets the last trajectory after its fork.  If that trajectory was the
  // only anomalous one, there are no anomalous trajectories after this call.
  void ResetLastSegment();

  // Deletes the last trajectory and removes it from |segments_|.  If there are
  // anomalous trajectories, their number is decremented and may become 0.
  void PopLastSegment();

  // Pops the burn of the man�uvre with the given index and all following
  // segments, then resets the last segment (which is the coast preceding
  // |man�uvres_[index]|).
  void PopSegmentsAffectedByMan�uvre(int index);

  // Reconstructs each man�uvre after |man�uvres_[index]| (starting with
  // |man�uvres_[index + 1]|), keeping the same burns but recomputing the
  // initial masses from |man�uvres_[index].final_mass()|.
  void UpdateInitialMassOfMan�uvresAfter(int index);

  Instant start_of_last_coast() const;

  // In the following functions, |index| refers to the index of a man�uvre.
  Instant start_of_burn(int index) const;
  Instant start_of_next_burn(int index) const;
  Instant start_of_previous_coast(int index) const;

  Mass const initial_mass_;
  Instant initial_time_;
  DegreesOfFreedom<Barycentric> initial_degrees_of_freedom_;
  Instant desired_final_time_;

  // The trajectory of the part, composed of any number of segments,
  // alternatively coasts and burns.
  DiscreteTrajectory<Barycentric> trajectory_;

  // Never empty; Starts and ends with a coast; coasts and burns alternate.
  std::vector<DiscreteTrajectorySegmentIterator<Barycentric>> segments_;
  // The last |anomalous_segments_| of |segments_| are anomalous, i.e., they
  // either end prematurely or follow an anomalous trajectory; in the latter
  // case they are empty.
  int anomalous_segments_ = 0;
  // The status of the first anomalous segment.  Set and used exclusively by
  // |ComputeSegments|.
  absl::Status anomalous_status_;

  std::vector<NavigationMan�uvre> man�uvres_;
  std::vector<not_null<std::unique_ptr<OrbitAnalyser>>> coast_analysers_;
  not_null<Ephemeris<Barycentric>*> ephemeris_;
  Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters_;
  Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
      generalized_adaptive_step_parameters_;
};

}  // namespace ksp_plugin
}  // namespace principia
