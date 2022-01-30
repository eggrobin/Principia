module;

#include "physics/rotating_body.hpp"
#include "quantities/named_quantities.hpp"

export module principia.ksp_plugin.equator_relevance_threshold;

import principia.ksp_plugin.frames;

namespace principia {

using quantities::Length;
using physics::RotatingBody;

export namespace ksp_plugin {

// Returns a distance from |body| that we consider is too far for the equator to
// be of interest.  Specifically, this distance is the maximum of
// - the semimajor axis of a supersynchronous orbit (1 orbit for 2 body
//   revolutions);
// - the distance at which a  |Geopotential| with a tolerance of 0x1p-24 starts
//   damping dynamical oblateness;
Length EquatorRelevanceThreshold(RotatingBody<Barycentric> const& body);

}  // namespace ksp_plugin
}  // namespace principia
