﻿#include "ksp_plugin/manœuvre_interface.hpp"

#include "quantities/constants.hpp"
#include "quantities/si.hpp"

namespace principia {

using constants::StandardGravity;
using geometry::RadiusLatitudeLongitude;
using si::Degree;
using si::Kilogram;
using si::Metre;
using si::Newton;
using si::Second;

namespace ksp_plugin {

namespace {

// TODO(egg): bad code replication! bad!

XYZ ToXYZ(R3Element<double> const& r3_element) {
  return {r3_element.x, r3_element.y, r3_element.z};
}

}  // namespace

Manœuvre* principia__NewManœuvreIspByWeight(
    double thrust,
    double initial_mass,
    double specific_impulse_by_weight,
    double right_ascension,
    double declination) {
  return std::make_unique<Manœuvre>(
      thrust * Newton,
      initial_mass * Kilogram,
      specific_impulse_by_weight * Second * StandardGravity,
      Vector<double, Frenet>(
          RadiusLatitudeLongitude(
              1.0,
              declination * Degree,
              right_ascension * Degree).ToCartesian())).release();
}

double principia__thrust(Manœuvre const* manœuvre) {
  return CHECK_NOTNULL(manœuvre)->thrust() / Newton;
}
double principia__initial_mass(Manœuvre const* manœuvre) {
  return CHECK_NOTNULL(manœuvre)->initial_mass() / Kilogram;
}
double principia__specific_impulse_by_weight(
    Manœuvre const* manœuvre) {
  return (CHECK_NOTNULL(manœuvre)->effective_exhaust_velocity() /
          StandardGravity) /
         Second;
}

XYZ principia__Δv(Manœuvre const* manœuvre) {
  CHECK_NOTNULL(manœuvre);
  return ToXYZ((manœuvre->direction() * manœuvre->Δv() /
                (Metre / Second)).coordinates());
}

void principia__set_duration(Manœuvre* manœuvre,
                             double duration) {
  CHECK_NOTNULL(manœuvre)->set_duration(duration * Second);
}

void principia__set_Δv(Manœuvre* manœuvre, double Δv) {
  CHECK_NOTNULL(manœuvre)->set_Δv(Δv * (Metre / Second));
}

double principia__initial_time(Manœuvre const* manœuvre) {
  return (CHECK_NOTNULL(manœuvre)->initial_time() - Instant()) / Second;
}

void principia__set_initial_time(Manœuvre* manœuvre,
                                 double initial_time) {
  return CHECK_NOTNULL(manœuvre)->
      set_initial_time(Instant(initial_time * Second));
}

double principia__time_of_half_Δv(Manœuvre const* manœuvre) {
  return (CHECK_NOTNULL(manœuvre)->time_of_half_Δv() - Instant()) / Second;
}

void principia__set_time_of_half_Δv(Manœuvre* manœuvre,
                                    double time_of_half_Δv) {
  return CHECK_NOTNULL(manœuvre)->
      set_time_of_half_Δv(Instant(time_of_half_Δv * Second));
}

}  // namespace ksp_plugin
}  // namespace principia
