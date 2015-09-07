#pragma once

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "base/not_null.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/frame_field.hpp"
#include "physics/massive_body.hpp"

namespace principia {

using base::not_null;

namespace physics {

// This class represent a pair of transformations of a trajectory from
// |FromFrame| to |ToFrame| with an intermediate representation in
// |ThroughFrame|.  Note that the trajectory in |ToFrame| is not the trajectory
// of a body since its past changes from moment to moment.
template<typename FromFrame, typename ThroughFrame, typename ToFrame>
class Transforms {
  static_assert(FromFrame::is_inertial && ToFrame::is_inertial,
                "Both FromFrame and ToFrame must be inertial");

 public:
  // A factory method where |ThroughFrame| is defined as follows: it has the
  // same axes as |FromFrame| and the body of |centre_trajectory| is the origin
  // of |ThroughFrame|.
  static not_null<std::unique_ptr<Transforms>> BodyCentredNonRotating(
      MassiveBody const& centre,
      ContinuousTrajectory<FromFrame> const& from_centre_trajectory,
      ContinuousTrajectory<ToFrame> const& to_centre_trajectory);

  // A factory method where |ThroughFrame| is defined as follows: its X axis
  // goes from the primary to the secondary bodies, its Y axis is in the plane
  // of the velocities of the bodies in their barycentric frame, on the same
  // side of the X axis as the velocity of the primary body, its Z axis is such
  // that it is right-handed.  The barycentre of the bodies is the origin of
  // |ThroughFrame|.
  static not_null<std::unique_ptr<Transforms>> BarycentricRotating(
      MassiveBody const& primary,
      ContinuousTrajectory<FromFrame> const& from_primary_trajectory,
      ContinuousTrajectory<ToFrame> const& to_primary_trajectory,
      MassiveBody const& secondary,
      ContinuousTrajectory<FromFrame> const& from_secondary_trajectory,
      ContinuousTrajectory<ToFrame> const& to_secondary_trajectory);

  // Use this only for testing!
  static not_null<std::unique_ptr<Transforms>> DummyForTesting();

  typename DiscreteTrajectory<FromFrame>::
  template TransformingIterator<ThroughFrame> first(
      DiscreteTrajectory<FromFrame> const& from_trajectory);

  typename DiscreteTrajectory<FromFrame>::
  template TransformingIterator<ThroughFrame> first_with_caching(
      not_null<DiscreteTrajectory<FromFrame>*> const from_trajectory);

  typename DiscreteTrajectory<FromFrame>::
  template TransformingIterator<ThroughFrame> first_on_or_after(
      DiscreteTrajectory<FromFrame> const& from_trajectory,
      Instant const& time);

  typename DiscreteTrajectory<ThroughFrame>::
  template TransformingIterator<ToFrame> second(
      Instant const& last,
      DiscreteTrajectory<ThroughFrame> const& through_trajectory);

  // The coordinate frame of |ThroughFrame|, expressed in the coordinates of
  // |ToFrame| at time |last|.
  FrameField<ToFrame> coordinate_frame(Instant const& last) const;

  // The acceleration on a test particle with the given |degrees_of_freedom|
  // at the given |time| due to the rotation and acceleration of |ThroughFrame|,
  // expressed in the coordinates of |ToFrame|.
  Vector<Acceleration, ToFrame> FrameAcceleration(
      DegreesOfFreedom<FromFrame> const& degrees_of_freedom,
      Instant const& time) const;

 private:
  // Just like a |DiscreteTrajectory::Transform|, except that the first
  // parameter is only bound when we know if we must cache.
  template<typename Frame1, typename Frame2>
  using FirstTransform = std::function<DegreesOfFreedom<Frame2>(
                             bool const,
                             Instant const&,
                             DegreesOfFreedom<Frame1> const&,
                             not_null<
                                 DiscreteTrajectory<Frame1> const*> const)>;
  FirstTransform<FromFrame, ThroughFrame> first_;

  // Just like a |DiscreteTrajectory::Transform|, except that the first
  // parameter is only bound when we know at what time (|last|) the transform
  // must be applied.
  template<typename Frame1, typename Frame2>
  using SecondTransform = std::function<DegreesOfFreedom<Frame2>(
                              Instant const& last,
                              Instant const& t,
                              DegreesOfFreedom<Frame1> const&,
                              not_null<
                                  DiscreteTrajectory<Frame1> const*> const)>;
  SecondTransform<ThroughFrame, ToFrame> second_;

  // A simple cache which monitors the hit rate.  Keyed by a pointer, so it
  // must get a notification when a trajectory is deleted.
  template<typename Frame1, typename Frame2>
  class Cache {
   public:
    bool Lookup(not_null<DiscreteTrajectory<Frame1> const*> const trajectory,
                Instant const& time,
                not_null<DegreesOfFreedom<Frame2>**> degrees_of_freedom);

    void Insert(not_null<DiscreteTrajectory<Frame1> const*> const trajectory,
                Instant const& time,
                DegreesOfFreedom<Frame2> const& degrees_of_freedom);

    void Delete(not_null<DiscreteTrajectory<Frame1> const*> const trajectory);

   private:
    std::map<not_null<DiscreteTrajectory<Frame1> const*>,
             std::map<Instant, DegreesOfFreedom<Frame2>>> cache_;
    std::map<not_null<DiscreteTrajectory<Frame1> const*>,
             std::int64_t> number_of_lookups_;
    std::map<not_null<DiscreteTrajectory<Frame1> const*>,
             std::int64_t> number_of_hits_;
  };

  // A cache for the result of the |first_| transform.  This cache assumes that
  // the iterator is never called with the same time but different degrees of
  // freedom.
  Cache<FromFrame, ThroughFrame> first_cache_;

  // Same as FrameField<ToFrame>, but the time is only bound when
  // |coordinate_frame| is called.
  std::function<Rotation<ToFrame, ToFrame>(
      Instant const& last, Position<ToFrame> const& q)> coordinate_frame_;

  std::function<Vector<Acceleration, ToFrame>(
      DegreesOfFreedom<FromFrame> const& degrees_of_freedom,
      Instant const& time)> frame_acceleration_;

  // Hints for the continuous trajectories.
  std::list<typename ContinuousTrajectory<FromFrame>::Hint> from_hints_;
  std::list<typename ContinuousTrajectory<ToFrame>::Hint> to_hints_;
};

}  // namespace physics
}  // namespace principia

#include "physics/transforms_body.hpp"
