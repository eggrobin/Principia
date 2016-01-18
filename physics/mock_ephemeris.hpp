
#pragma once

#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock-generated-function-mockers.h"
#include "gmock/gmock-matchers.h"
#include "gmock/gmock.h"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
template <typename DifferentialEquation> class AdaptiveStepSizeIntegrator;
template <typename DifferentialEquation> class FixedStepSizeIntegrator;
}  // namespace integrators
namespace serialization {
class Ephemeris;
}  // namespace serialization
}  // namespace principia

namespace principia {
namespace physics {

class MassiveBody;
template <typename Frame> class ContinuousTrajectory;

template<typename Frame>
class MockEphemeris : public Ephemeris<Frame> {
 public:
  MockEphemeris() : Ephemeris<Frame>() {}

  MOCK_CONST_METHOD0_T(bodies,
                       std::vector<not_null<MassiveBody const*>> const&());
  MOCK_CONST_METHOD1_T(trajectory,
                       not_null<ContinuousTrajectory<Frame> const*>(
                           not_null<MassiveBody const*> body));
  MOCK_CONST_METHOD0_T(empty, bool());
  MOCK_CONST_METHOD0_T(t_min, Instant());
  MOCK_CONST_METHOD0_T(t_max, Instant());
  MOCK_CONST_METHOD0_T(
      planetary_integrator,
      FixedStepSizeIntegrator<
          typename Ephemeris<Frame>::NewtonianMotionEquation> const&());

  MOCK_METHOD1_T(ForgetBefore, void(Instant const& t));
  MOCK_METHOD1_T(Prolong, void(Instant const& t));
  MOCK_METHOD6_T(
      FlowWithAdaptiveStep,
      void(not_null<DiscreteTrajectory<Frame>*> const trajectory,
           typename Ephemeris<Frame>::IntrinsicAcceleration
               intrinsic_acceleration,
           Length const& length_integration_tolerance,
           Speed const& speed_integration_tolerance,
           AdaptiveStepSizeIntegrator<
               typename Ephemeris<Frame>::NewtonianMotionEquation> const&
               integrator,
           Instant const& t));
  MOCK_METHOD4_T(
      FlowWithFixedStep,
      void(std::vector<not_null<DiscreteTrajectory<Frame>*>> const&
               trajectories,
           typename Ephemeris<Frame>::IntrinsicAccelerations const&
               intrinsic_accelerations,
           Time const& step,
           Instant const& t));

  MOCK_CONST_METHOD2_T(
      ComputeGravitationalAccelerationOnMasslessBody,
      Vector<Acceleration, Frame>(Position<Frame> const& position,
                                  Instant const & t));

  // NOTE(phl): This overload introduces ambiguities in the expectations.
  // MOCK_CONST_METHOD2_T(
  //     ComputeGravitationalAccelerationOnMasslessBody,
  //     Vector<Acceleration, Frame>(
  //         not_null<DiscreteTrajectory<Frame>*> /*const*/ trajectory,
  //         Instant const& t));

  // NOTE(phl): The commented-out const below is to work-around a compiler
  // internal error.  Don't ask.
  MOCK_CONST_METHOD2_T(
      ComputeGravitationalAccelerationOnMassiveBody,
      Vector<Acceleration, Frame>(
          not_null<MassiveBody const*> /*const*/ body,
          Instant const& t));

  MOCK_CONST_METHOD1_T(serialization_index_for_body,
                       int(not_null<MassiveBody const*> const body));
  MOCK_CONST_METHOD1_T(
      body_for_serialization_index,
      not_null<MassiveBody const*>(int const serialization_index));

  MOCK_CONST_METHOD1_T(WriteToMessage,
                       void(not_null<serialization::Ephemeris*> const message));
};

}  // namespace physics
}  // namespace principia
