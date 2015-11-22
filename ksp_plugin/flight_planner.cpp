#include "ksp_plugin/flight_planner.hpp"

#include "ksp_plugin/interface.hpp"
#include "quantities/constants.hpp"

namespace principia {

using quantities::constants::StandardGravity;
using quantities::si::Kilo;
using quantities::si::Newton;

namespace ksp_plugin {

Burn::Burn(BurnInterface burn_interface)
    : name(burn_interface.name),
      specific_impulse(
          burn_interface.specific_impulse * (Second * StandardGravity)),
      thrust(burn_interface.thrust * Kilo(Newton)),
      Δv(Velocity<Frenet<Rendering>>(ToR3Element(burn_interface.Δv) *
             (Metre / Second))),
      initial_time_from_end_of_previous(
          burn_interface.initial_time_from_end_of_previous * Second) {}

}  // namespace ksp_plugin
}  // namespace principia
