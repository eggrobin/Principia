module;

#include "physics/discrete_trajectory_segment.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

export module principia.ksp_plugin.integrators;

import principia.ksp_plugin.frames;

namespace principia {

using physics::DiscreteTrajectorySegment;
using physics::Ephemeris;
using quantities::Length;
using quantities::si::Metre;
using quantities::si::Milli;

export namespace ksp_plugin {

// Parameters for downsampling after fixed-step integration.
DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters
DefaultDownsamplingParameters();

// Factories for parameters used to control integration.
Ephemeris<Barycentric>::AccuracyParameters
DefaultEphemerisAccuracyParameters();
Ephemeris<Barycentric>::FixedStepParameters
DefaultEphemerisFixedStepParameters();
Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
DefaultBurnParameters();
Ephemeris<Barycentric>::FixedStepParameters DefaultHistoryParameters();
Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPredictionParameters();
Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPsychohistoryParameters();

}  // namespace ksp_plugin
}  // namespace principia
