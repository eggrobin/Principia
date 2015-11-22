#pragma once

#include <deque>
#include <string>

#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {
namespace ksp_plugin {

struct BurnInterface {
  char const* const name;
  // In s * g_0.
  double const specific_impulse;
  // In kN.
  double const thrust;
  // In m/s.
  XYZ const ?v;
  // In s.
  double const initial_time_from_end_of_previous;
};

static_assert(std::is_standard_layout<BurnInterface>::value,
              "BurnInterface is used for interfacing");

struct Burn {
  std::string const name;
  SpecificImpulse const specific_impulse;
  Force const thrust;
  Velocity<Frenet<Rendering>> const ?v;
  Time const initial_time_from_end_of_previous;

  explicit Burn(BurnInterface burn_interface);
};

struct FlightPlanner {
  Mass const initial_mass;
  Instant const start;
  Instant end;
  std::deque<Burn> burns;
};

extern "C" DLLEXPORT
int CDECL AppendBurn(FlightPlanner* const planner, BurnInterface const burn);
extern "C" DLLEXPORT
void CDECL DeleteLastBurn(FlightPlanner* const planner, int const index);
extern "C" DLLEXPORT
void CDECL EditBurn(FlightPlanner* const planner,
                    int const index,
                    BurnInterface const burn);

}  // namespace ksp_plugin
}  // namespace principia
