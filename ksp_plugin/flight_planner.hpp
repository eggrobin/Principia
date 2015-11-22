#pragma once

#include <deque>
#include <string>

#include "ksp_plugin/frames.hpp"
#include "physics/dynamic_frame.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using physics::DynamicFrame;
using physics::Frenet;
using quantities::Force;
using quantities::Mass;
using quantities::SpecificImpulse;
using quantities::Time;

namespace ksp_plugin {

struct BurnInterface;

struct Burn {
  std::string name;
  SpecificImpulse specific_impulse;
  Force thrust;
  Velocity<Frenet<Rendering>> Δv;
  Time initial_time_from_end_of_previous;

  explicit Burn(BurnInterface burn_interface);
};

struct FlightPlanner {
  Mass const initial_mass;
  Instant const initial_time;
  Instant final_time;
  std::deque<Burn> burns;
  DynamicFrame<Barycentric, Rendering> const* frame;
};

}  // namespace ksp_plugin
}  // namespace principia
