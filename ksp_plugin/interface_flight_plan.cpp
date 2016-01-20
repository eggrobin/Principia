﻿
#include <memory>
#include <ostream>
#include <type_traits>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/epoch.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "journal/method.hpp"
#include "journal/profiles.hpp"
#include "ksp_plugin/burn.hpp"
#include "ksp_plugin/flight_plan.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "ksp_plugin/plugin.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/massive_body.hpp"
#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::not_null;
using geometry::Instant;
using ksp_plugin::Barycentric;
using ksp_plugin::FlightPlan;
using ksp_plugin::Navigation;
using ksp_plugin::NavigationManœuvre;
using ksp_plugin::Vessel;
using quantities::constants::StandardGravity;
using quantities::si::Kilo;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Newton;
using quantities::si::Second;
using quantities::si::Tonne;

namespace interface {

namespace {

not_null<Vessel*> GetVessel(Plugin const* const plugin,
                            char const* const vessel_guid) {
  // TODO(phl): Should this check if the vessel is synchronized/initialized?
  CHECK(CHECK_NOTNULL(plugin)->HasVessel(vessel_guid)) << vessel_guid;
  return plugin->GetVessel(vessel_guid);
}

FlightPlan& GetFlightPlan(Plugin const* const plugin,
                          char const* const vessel_guid) {
  Vessel const& vessel = *GetVessel(plugin, vessel_guid);
  CHECK(vessel.has_flight_plan()) << vessel_guid;
  return *vessel.flight_plan();
}

Burn GetBurn(NavigationManœuvre const& manœuvre) {
  return {manœuvre.thrust() / Kilo(Newton),
          manœuvre.specific_impulse() / (Second * StandardGravity),
          principia__GetNavigationFrameParameters(manœuvre.frame()),
          (manœuvre.initial_time() - Instant()) / Second,
          ToXYZ((manœuvre.Δv() / (Metre / Second)) *
                    manœuvre.direction().coordinates())};
}

ksp_plugin::Burn ToBurn(Plugin const* const plugin, Burn const& burn) {
  return {burn.thrust_in_kilonewtons * Kilo(Newton),
          burn.specific_impulse_in_seconds_g0 * Second * StandardGravity,
          base::check_not_null(std::unique_ptr<NavigationFrame>(
              principia__NewNavigationFrame(plugin, burn.frame))),
          Instant() + burn.initial_time * Second,
          Velocity<Frenet<Navigation>>(
              ToR3Element(burn.delta_v) * (Metre / Second))};
}

NavigationManoeuvre ToNavigationManoeuvre(NavigationManœuvre const& manœuvre) {
  NavigationManoeuvre result;
  result.burn = GetBurn(manœuvre);
  result.initial_mass_in_tonnes = manœuvre.initial_mass() / Tonne;
  result.final_mass_in_tonnes = manœuvre.final_mass() / Tonne;
  result.mass_flow = manœuvre.mass_flow() / (Kilogram / Second);
  result.duration = manœuvre.duration() / Second;
  result.final_time = (manœuvre.final_time() - Instant()) / Second;
  result.time_of_half_delta_v =
      (manœuvre.time_of_half_Δv() - Instant()) / Second;
  result.time_to_half_delta_v = manœuvre.time_to_half_Δv() / Second;
  result.direction = ToXYZ(manœuvre.direction().coordinates());
  return result;
}

}  // namespace

bool principia__FlightPlanAppend(Plugin const* const plugin,
                                 char const* const vessel_guid,
                                 Burn const burn) {
  journal::Method<journal::FlightPlanAppend> m({plugin, vessel_guid, burn});
  return m.Return(GetFlightPlan(plugin, vessel_guid).
                      Append(ToBurn(plugin, burn)));
}

void principia__FlightPlanCreate(Plugin const* const plugin,
                                 char const* const vessel_guid,
                                 double const final_time,
                                 double const mass_in_tonnes) {
  journal::Method<journal::FlightPlanCreate> m({plugin,
                                                vessel_guid,
                                                final_time,
                                                mass_in_tonnes});
  CHECK_NOTNULL(plugin)->CreateFlightPlan(vessel_guid,
                                          Instant() + final_time * Second,
                                          mass_in_tonnes * Tonne);
  return m.Return();
}

void principia__FlightPlanDelete(Plugin const* const plugin,
                                 char const* const vessel_guid) {
  journal::Method<journal::FlightPlanDelete> m({plugin, vessel_guid});
  GetVessel(plugin, vessel_guid)->DeleteFlightPlan();
  return m.Return();
}

NavigationManoeuvre principia__FlightPlanGetManoeuvre(
    Plugin const* const plugin,
    char const* const vessel_guid,
    int const index) {
  journal::Method<journal::FlightPlanGetManoeuvre> m({plugin,
                                                      vessel_guid,
                                                      index});
  return m.Return(ToNavigationManoeuvre(
             GetFlightPlan(plugin, vessel_guid).GetManœuvre(index)));
}

int principia__FlightPlanNumberOfManoeuvres(Plugin const* const plugin,
                                            char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfManoeuvres> m({plugin,
                                                            vessel_guid});
  return m.Return(GetFlightPlan(plugin, vessel_guid).number_of_manœuvres());
}

int principia__FlightPlanNumberOfSegments(Plugin const* const plugin,
                                          char const* const vessel_guid) {
  journal::Method<journal::FlightPlanNumberOfSegments> m({plugin, vessel_guid});
  return m.Return(GetFlightPlan(plugin, vessel_guid).number_of_segments());
}

void principia__FlightPlanRemoveLast(Plugin const* const plugin,
                                     char const* const vessel_guid) {
  journal::Method<journal::FlightPlanRemoveLast> m({plugin, vessel_guid});
  GetFlightPlan(plugin, vessel_guid).RemoveLast();
  return m.Return();
}

LineAndIterator* principia__FlightPlanRenderedSegment(
    Plugin const* const plugin,
    char const* const vessel_guid,
    XYZ const sun_world_position,
    int const index) {
  journal::Method<journal::FlightPlanRenderedSegment> m({plugin,
                                                         vessel_guid,
                                                         sun_world_position,
                                                         index});
  DiscreteTrajectory<Barycentric>::Iterator begin;
  DiscreteTrajectory<Barycentric>::Iterator end;
  GetFlightPlan(plugin, vessel_guid).GetSegment(index, &begin, &end);
  RenderedTrajectory<World> rendered_trajectory = CHECK_NOTNULL(plugin)->
      RenderedTrajectoryFromIterators(
          begin, end,
          World::origin + Displacement<World>(
                              ToR3Element(sun_world_position) * Metre));
  not_null<std::unique_ptr<LineAndIterator>> result =
      make_not_null_unique<LineAndIterator>(std::move(rendered_trajectory));
  result->it = result->rendered_trajectory.begin();
  return m.Return(result.release());
}

bool principia__FlightPlanReplaceLast(Plugin const* const plugin,
                                      char const* const vessel_guid,
                                      Burn const burn) {
  journal::Method<journal::FlightPlanReplaceLast> m({plugin,
                                                     vessel_guid,
                                                     burn});
  return m.Return(GetFlightPlan(plugin, vessel_guid).
                      ReplaceLast(ToBurn(plugin, burn)));
}

bool principia__FlightPlanSetFinalTime(Plugin const* const plugin,
                                       char const* const vessel_guid,
                                       double const final_time) {
  journal::Method<journal::FlightPlanSetFinalTime> m({plugin,
                                                      vessel_guid,
                                                      final_time});
  return m.Return(GetFlightPlan(plugin, vessel_guid).
                      SetFinalTime(Instant() + final_time * Second));
}

void principia__FlightPlanSetTolerances(
    Plugin const* const plugin,
    char const* const vessel_guid,
    double const length_integration_tolerance,
    double const speed_integration_tolerance) {
  journal::Method<journal::FlightPlanSetTolerances>
      m({plugin,
         vessel_guid,
         length_integration_tolerance,
         speed_integration_tolerance});
  GetFlightPlan(plugin, vessel_guid).
      SetTolerances(length_integration_tolerance * Metre,
                    speed_integration_tolerance * (Metre / Second));
  return m.Return();
}

}  // namespace interface
}  // namespace principia
