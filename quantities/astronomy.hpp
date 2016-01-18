﻿
#pragma once

#include "geometry/epoch.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massive_body.hpp"
#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

// This namespace contains units commonly used in astronomy that are not
// accepted for use with the SI.
namespace astronomy {

quantities::Mass   const SolarMass     = 1.98855e30 * si::Kilogram;
quantities::Mass   const JupiterMass   = 1.8986e27 * si::Kilogram;
quantities::Mass   const EarthMass     = 5.9742e24 * si::Kilogram;
quantities::Time   const JulianYear    = 365.25 * si::Day;
quantities::Length const Parsec        = 648000 / π * si::AstronomicalUnit;
quantities::Length const LightYear     = constants::SpeedOfLight * JulianYear;
quantities::Length const LunarDistance = 384400000 * si::Metre;

}  // namespace astronomy
}  // namespace quantities
}  // namespace principia
