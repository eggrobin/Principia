
#pragma once

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock-generated-function-mockers.h"
#include "gmock/gmock-matchers.h"
#include "gmock/gmock.h"
#include "ksp_plugin/burn.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {

class MockFlightPlan : public FlightPlan {
 public:
  MockFlightPlan() : FlightPlan() {}

  MOCK_CONST_METHOD0(number_of_manœuvres, int());
  MOCK_CONST_METHOD1(GetManœuvre, NavigationManœuvre const& (int const index));

  MOCK_METHOD0(RemoveLast, void());

  MOCK_CONST_METHOD1(AppendConstRef, bool(Burn const& burn));
  MOCK_CONST_METHOD1(ReplaceLastConstRef, bool(Burn const& burn));

  bool Append(Burn burn);
  bool ReplaceLast(Burn burn);

  MOCK_METHOD1(SetFinalTime, bool(Instant const& final_time));

  MOCK_METHOD2(SetTolerances,
               void(Length const& length_integration_tolerance,
                    Speed const& speed_integration_tolerance));

  MOCK_CONST_METHOD0(number_of_segments, int());

  MOCK_CONST_METHOD3(
      GetSegment,
      void(int const index,
           not_null<DiscreteTrajectory<Barycentric>::Iterator*> begin,
           not_null<DiscreteTrajectory<Barycentric>::Iterator*> end));
};

}  // namespace ksp_plugin
}  // namespace principia
