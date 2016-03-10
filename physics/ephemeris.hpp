﻿
#pragma once

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "google/protobuf/repeated_field.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {

using geometry::Position;
using geometry::Vector;
using integrators::AdaptiveStepSizeIntegrator;
using integrators::FixedStepSizeIntegrator;
using integrators::SpecialSecondOrderDifferentialEquation;

namespace physics {

template<typename Frame>
class Ephemeris {
  static_assert(Frame::is_inertial, "Frame must be inertial");

 public:
  using IntrinsicAcceleration =
      std::function<Vector<Acceleration, Frame>(Instant const& time)>;
  static std::nullptr_t constexpr kNoIntrinsicAcceleration = nullptr;
  using IntrinsicAccelerations = std::vector<IntrinsicAcceleration>;
  static IntrinsicAccelerations const kNoIntrinsicAccelerations;

  // The equation describing the motion of the |bodies_|.
  using NewtonianMotionEquation =
      SpecialSecondOrderDifferentialEquation<Position<Frame>>;

  class AdaptiveStepParameters {
   public:
    // The |length_| and |speed_integration_tolerance|s are used to compute the
    // |tolerance_to_error_ratio| for step size control.  The number of steps is
    // limited to |max_steps|.
    AdaptiveStepParameters(
        AdaptiveStepSizeIntegrator<NewtonianMotionEquation> const& integrator,
        std::int64_t const max_steps,
        Length const& length_integration_tolerance,
        Speed const& speed_integration_tolerance);

    // These functions can serialize to/from any message having the right
    // fields.  It would be nicer if this was a separate message or a part of
    // the Ephemeris message, but that would break pre-Буняко́вский
    // compatibility.
    // TODO(phl): Restructure this after Cantor.
    template<typename T>
    void WriteToMessage(not_null<T*> const t) const;
    template<typename T>
    static AdaptiveStepParameters ReadFromMessage(T const& t);

   private:
    // This will refer to a static object returned by a factory.
    not_null<AdaptiveStepSizeIntegrator<NewtonianMotionEquation> const*>
        integrator_;
    std::int64_t max_steps_;
    Length length_integration_tolerance_;
    Speed speed_integration_tolerance_;
    friend class Ephemeris<Frame>;
  };

  class FixedStepParameters {
   public:
    FixedStepParameters(
        FixedStepSizeIntegrator<NewtonianMotionEquation> const& integrator,
        Time const& step);

    Time const& step() const;

   private:
    // This will refer to a static object returned by a factory.
    not_null<FixedStepSizeIntegrator<NewtonianMotionEquation> const*>
        integrator_;
    Time step_;
    friend class Ephemeris<Frame>;
  };

  // Constructs an Ephemeris that owns the |bodies|.  The elements of vectors
  // |bodies| and |initial_state| correspond to one another.
  Ephemeris(std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies,
            std::vector<DegreesOfFreedom<Frame>> const& initial_state,
            Instant const& initial_time,
            Length const& fitting_tolerance,
            FixedStepParameters const& parameters);

  virtual ~Ephemeris() = default;

  // Returns the bodies in the order in which they were given at construction.
  virtual std::vector<not_null<MassiveBody const*>> const& bodies() const;

  // Returns the trajectory for the given |body|.
  virtual not_null<ContinuousTrajectory<Frame> const*> trajectory(
      not_null<MassiveBody const*> const body) const;

  // Returns true if at least one of the trajectories is empty.
  virtual bool empty() const;

  // The maximum of the |t_min|s of the trajectories.
  virtual Instant t_min() const;
  // The mimimum of the |t_max|s of the trajectories.
  virtual Instant t_max() const;

  virtual FixedStepSizeIntegrator<NewtonianMotionEquation> const&
  planetary_integrator() const;

  // Calls |ForgetAfter| on all trajectories for a time which is greater than or
  // equal to |t|, and less than 6 months after |t|.  On return |t_max() >= t|.
  virtual void ForgetAfter(Instant const& t);

  // Calls |ForgetBefore| on all trajectories.  On return |t_min() == t|.
  virtual void ForgetBefore(Instant const& t);

  // Prolongs the ephemeris up to at least |t|.  After the call, |t_max() >= t|.
  virtual void Prolong(Instant const& t);

  // Integrates, until exactly |t| (except for timeouts or singularities), the
  // |trajectory| followed by a massless body in the gravitational potential
  // described by |*this|.  If |t > t_max()|, calls |Prolong(t)| beforehand.
  // Returns true if and only if |*trajectory| was integrated until |t|.
  virtual bool FlowWithAdaptiveStep(
      not_null<DiscreteTrajectory<Frame>*> const trajectory,
      IntrinsicAcceleration intrinsic_acceleration,
      Instant const& t,
      AdaptiveStepParameters const& parameters);

  // Integrates, until at most |t|, the |trajectories| followed by massless
  // bodies in the gravitational potential described by |*this|.  If
  // |t > t_max()|, calls |Prolong(t)| beforehand.
  virtual void FlowWithFixedStep(
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
      IntrinsicAccelerations const& intrinsic_accelerations,
      Instant const& t,
      FixedStepParameters const& parameters);

  // Returns the gravitational acceleration on a massless body located at the
  // given |position| at time |t|.
  virtual Vector<Acceleration, Frame>
  ComputeGravitationalAccelerationOnMasslessBody(
      Position<Frame> const& position,
      Instant const & t) const;

  // Returns the gravitational acceleration on the massless body having the
  // given |trajectory| at time |t|.  |t| must be one of the times of the
  // |trajectory|.
  virtual Vector<Acceleration, Frame>
  ComputeGravitationalAccelerationOnMasslessBody(
      not_null<DiscreteTrajectory<Frame>*> const trajectory,
      Instant const& t) const;

  // Returns the gravitational acceleration on the massive |body| at time |t|.
  // |body| must be one of the bodies of this object.
  virtual Vector<Acceleration, Frame>
  ComputeGravitationalAccelerationOnMassiveBody(
      not_null<MassiveBody const*> const body,
      Instant const& t) const;

  // Returns the index of the given body in the serialization produced by
  // |WriteToMessage| and read by the |Read...| functions.  This index is not
  // suitable for other uses.
  virtual int serialization_index_for_body(
      not_null<MassiveBody const*> const body) const;

  virtual not_null<MassiveBody const*> body_for_serialization_index(
      int const serialization_index) const;

  virtual void WriteToMessage(
      not_null<serialization::Ephemeris*> const message) const;
  static not_null<std::unique_ptr<Ephemeris>> ReadFromMessage(
      serialization::Ephemeris const& message);

  // Compatibility method for construction an ephemeris from pre-Bourbaki data.
  static std::unique_ptr<Ephemeris> ReadFromPreBourbakiMessages(
      google::protobuf::RepeatedPtrField<
          serialization::Plugin::CelestialAndProperties> const& messages,
      Length const& fitting_tolerance,
      typename Ephemeris<Frame>::FixedStepParameters const& fixed_parameters);

 protected:
  // For mocking purposes, leaves everything uninitialized and uses a dummy
  // integrator that crashes when used.
  Ephemeris();

 private:
  void AppendMassiveBodiesState(
      typename NewtonianMotionEquation::SystemState const& state);
  static void AppendMasslessBodiesState(
      typename NewtonianMotionEquation::SystemState const& state,
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories);

  // Computes the acceleration due to one body, |body1| (with index |b1| in the
  // |positions| and |accelerations| arrays) on the bodies |bodies2| (with
  // indices [b2_begin, b2_end[ in the |positions| and |accelerations| arrays).
  // The template parameters specify what we know about the bodies, and
  // therefore what forces apply.
  template<bool body1_is_oblate,
           bool body2_is_oblate>
  static void ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies(
      MassiveBody const& body1,
      size_t const b1,
      std::vector<not_null<MassiveBody const*>> const& bodies2,
      size_t const b2_begin,
      size_t const b2_end,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations);

  // Computes the accelerations due to one body, |body1| (with index |b1| in the
  // |hints|, |bodies_| and |trajectories_| arrays) on massless bodies at the
  // given |positions|.  The template parameter specifies what we know about the
  // massive body, and therefore what forces apply.
  template<bool body1_is_oblate>
  void ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies(
      Instant const& t,
      MassiveBody const& body1,
      size_t const b1,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations,
      not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
          const hints) const;

  // Computes the accelerations between all the massive bodies in |bodies_|.
  void ComputeMassiveBodiesGravitationalAccelerations(
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const
          accelerations) const;

  // Computes the acceleration exerted by the massive bodies in |bodies_| on
  // massless bodies.  The massless bodies are at the given |positions|.  The
  // |hints| are passed to
  // ComputeGravitationalAccelerationByMassiveBodyOnMasslessBody for efficient
  // computation of the positions of the massive bodies.
  void ComputeMasslessBodiesGravitationalAccelerations(
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations,
      not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
          const hints) const;

  // Same as above, but the massless bodies have intrinsic accelerations.
  // |intrinsic_accelerations| may be empty.
  void ComputeMasslessBodiesTotalAccelerations(
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
      std::vector<IntrinsicAcceleration> const& intrinsic_accelerations,
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations,
      not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
          const hints) const;

  // Computes an estimate of the ratio |tolerance / error|.
  static double ToleranceToErrorRatio(
      Length const& length_integration_tolerance,
      Speed const& speed_integration_tolerance,
      Time const& current_step_size,
      typename NewtonianMotionEquation::SystemStateError const& error);

  // The bodies in the order in which they were given at construction.
  std::vector<not_null<MassiveBody const*>> unowned_bodies_;

  // The indices of bodies in |unowned_bodies_|.
  std::map<not_null<MassiveBody const*>, int> unowned_bodies_indices_;

  // The oblate bodies precede the spherical bodies in this vector.  The system
  // state is indexed in the same order.
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies_;

  // The indices in |bodies_| correspond to those in |oblate_bodies_| and
  // |spherical_bodies_|, in sequence.  The elements of |oblate_bodies_| are
  // really |OblateBody<Frame>| but it's inconvenient to express.
  std::vector<not_null<MassiveBody const*>> oblate_bodies_;
  std::vector<not_null<MassiveBody const*>> spherical_bodies_;

  // The indices in |bodies_| correspond to those in |trajectories_|.
  std::vector<not_null<ContinuousTrajectory<Frame>*>> trajectories_;

  std::map<not_null<MassiveBody const*>,
           not_null<std::unique_ptr<ContinuousTrajectory<Frame>>>>
      bodies_to_trajectories_;

  FixedStepParameters const parameters_;
  Length const fitting_tolerance_;
  typename NewtonianMotionEquation::SystemState last_state_;

  // These are the states other that the last which we preserve in order to be
  // to implement ForgetAfter.  The |state.time.value| are |t_max()| values for
  // all the underlying trajectories.
  std::vector<typename NewtonianMotionEquation::SystemState>
      intermediate_states_;

  int number_of_oblate_bodies_ = 0;
  int number_of_spherical_bodies_ = 0;

  NewtonianMotionEquation massive_bodies_equation_;
};

}  // namespace physics
}  // namespace principia

#include "physics/ephemeris_body.hpp"
