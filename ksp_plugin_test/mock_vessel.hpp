﻿
#pragma once

#include "gmock/gmock.h"
#include "ksp_plugin/vessel.hpp"

namespace principia {
namespace ksp_plugin {

class MockVessel : public Vessel {
 public:
  MockVessel() : Vessel() {}

  MOCK_CONST_METHOD0(body, not_null<MasslessBody const*>());
  MOCK_CONST_METHOD0(is_initialized, bool());

  MOCK_CONST_METHOD0(parent, not_null<Celestial const*>());
  MOCK_METHOD1(set_parent, void(not_null<Celestial const*> const parent));

  MOCK_CONST_METHOD0(prolongation, DiscreteTrajectory<Barycentric> const&());
  MOCK_METHOD0(mutable_prolongation,
               not_null<DiscreteTrajectory<Barycentric>*>());

  MOCK_CONST_METHOD0(flight_plan, not_null<FlightPlan*>());
  MOCK_CONST_METHOD0(has_flight_plan, bool());

  MOCK_CONST_METHOD0(prediction, DiscreteTrajectory<Barycentric> const&());
  MOCK_CONST_METHOD0(has_prediction, bool());

  MOCK_METHOD0(set_dirty, void());
  MOCK_CONST_METHOD0(is_dirty, bool());

  MOCK_METHOD2(CreateHistoryAndForkProlongation,
               void(Instant const& time,
                    DegreesOfFreedom<Barycentric> const& degrees_of_freedom));

  MOCK_METHOD1(AdvanceTime, void(Instant const& time));

  MOCK_METHOD1(ForgetBefore, void(Instant const& time));

  MOCK_CONST_METHOD0(ForgettableTime, Instant());

  MOCK_METHOD3(CreateFlightPlan,
               void(Instant const& final_time,
                    Mass const& initial_mass,
                    Ephemeris<Barycentric>::AdaptiveStepParameters const&
                        adaptive_parameters));

  MOCK_METHOD0(DeleteFlightPlan, void());

  MOCK_METHOD2(UpdatePrediction,
               void(Instant const& last_time,
                    Ephemeris<Barycentric>::AdaptiveStepParameters const&
                        adaptive_parameters));

  MOCK_METHOD0(DeletePrediction, void());

  MOCK_CONST_METHOD1(WriteToMessage, void(
      not_null<serialization::Vessel*> const message));
};

}  // namespace ksp_plugin
}  // namespace principia
