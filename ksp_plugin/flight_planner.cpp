#include "ksp_plugin/flight_planner.hpp"

#include "quantities/constants.hpp"

namespace principia {

using quantities::constants::StandardGravity;
using quantities::si::Kilo;
using quantities::si::Newton;

namespace ksp_plugin {

namespace {

R3Element<double> ToR3Element(XYZ const& xyz) {
  return {xyz.x, xyz.y, xyz.z};
}

}  // namespace

Burn::Burn(BurnInterface burn_interface)
    : name(burn_interface.name),
      specific_impulse(
          burn_interface.specific_impulse * (Second * StandardGravity)),
      thrust(burn_interface.thrust * Kilo(Newton)),
      Δv(Velocity<Frenet<Rendering>>(ToR3Element(burn_interface.Δv) *
             (Metre / Second))),
      initial_time_from_end_of_previous(
          burn_interface.initial_time_from_end_of_previous * Second) {}

int AppendBurn(FlightPlanner* const planner, BurnInterface const burn) {
  planner->burns.emplace_back(Burn(burn));
  return planner->burns.size() - 1;
}

void DeleteLastBurn(FlightPlanner* const planner, int const index) {
  planner->burns.pop_back();
  CHECK_EQ(planner->burns.size(), index);
}

void EditBurn(FlightPlanner* const planner,
              int const index,
              BurnInterface const burn) {
  planner->burns[index] = Burn(burn);
}

}  // namespace ksp_plugin
}  // namespace principia
