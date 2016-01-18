#pragma once

#include <memory>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/forkable.hpp"
#include "physics/massive_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {

namespace geometry {
template <typename FrameTag, FrameTag frame_tag, bool frame_is_inertial> class Frame;
}  // namespace geometry

using geometry::Instant;
using geometry::Velocity;
using physics::Frenet;
using quantities::Force;
using quantities::SpecificImpulse;

namespace ksp_plugin {

// Parameters for constructing a |NavigationManœuvre|, excluding the initial
// mass.  This owns a |NavigationFrame| and is therefore not copyable.
struct Burn {
  Force thrust;
  SpecificImpulse specific_impulse;
  not_null<std::unique_ptr<NavigationFrame const>> frame;
  Instant initial_time;
  Velocity<Frenet<Navigation>> Δv;
};

NavigationManœuvre MakeNavigationManœuvre(Burn burn, Mass const& initial_mass);

}  // namespace ksp_plugin
}  // namespace principia
