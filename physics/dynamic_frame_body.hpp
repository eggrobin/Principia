
#pragma once

#include <memory>
#include <ostream>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/epoch.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/rotation.hpp"
#include "glog/logging.h"
#include "google/protobuf/extension_set.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/barycentric_rotating_dynamic_frame.hpp"
#include "physics/body_centered_non_rotating_dynamic_frame.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/physics.pb.h"

namespace principia {

using geometry::Bivector;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::R3x3Matrix;
using geometry::Wedge;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Second;

namespace physics {

template<typename InertialFrame, typename ThisFrame>
Rotation<Frenet<ThisFrame>, ThisFrame>
DynamicFrame<InertialFrame, ThisFrame>::FrenetFrame(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  Velocity<ThisFrame> const& velocity = degrees_of_freedom.velocity();
  Vector<Acceleration, ThisFrame> const acceleration =
      GeometricAcceleration(t, degrees_of_freedom);
  Vector<Acceleration, ThisFrame> normal_acceleration = acceleration;
  velocity.template Orthogonalize<Acceleration>(&normal_acceleration);
  Vector<double, ThisFrame> tangent = Normalize(velocity);
  Vector<double, ThisFrame> normal = Normalize(normal_acceleration);
  Bivector<double, ThisFrame> binormal = Wedge(tangent, normal);
  // Maps |tangent| to {1, 0, 0}, |normal| to {0, 1, 0}, and |binormal| to
  // {0, 0, 1}.
  return Rotation<Frenet<ThisFrame>, ThisFrame>(
      R3x3Matrix(tangent.coordinates(),
                 normal.coordinates(),
                 binormal.coordinates()).Transpose());
}

template<typename InertialFrame, typename ThisFrame>
std::unique_ptr<DynamicFrame<InertialFrame, ThisFrame>>
DynamicFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::DynamicFrame const& message) {
  std::unique_ptr<DynamicFrame> result;
  int extensions_found = 0;
  // TODO(egg): See if the reset/release pairs could be avoided by smart
  // conversions of |not_null|.
  if (message.HasExtension(serialization::BarycentricRotatingDynamicFrame::
                               barycentric_rotating_dynamic_frame)) {
    ++extensions_found;
    result.reset(BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
        ReadFromMessage(
            ephemeris,
            message.GetExtension(
              serialization::BarycentricRotatingDynamicFrame::
                  barycentric_rotating_dynamic_frame)).release());
  }
  if (message.HasExtension(serialization::BodyCentredNonRotatingDynamicFrame::
                               body_centred_non_rotating_dynamic_frame)) {
    ++extensions_found;
    result.reset(BodyCentredNonRotatingDynamicFrame<InertialFrame, ThisFrame>::
        ReadFromMessage(
            ephemeris,
            message.GetExtension(
                serialization::BodyCentredNonRotatingDynamicFrame::
                    body_centred_non_rotating_dynamic_frame)).release());
  }
  CHECK_LE(extensions_found, 1) << message.DebugString();
  // For pre-Brouwer compatibility, return a null pointer if no extension is
  // found.
  return result;
}

}  // namespace physics
}  // namespace principia
