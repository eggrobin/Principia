
#pragma once

#include <ostream>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "ksp_plugin/celestial.hpp"
#include "mathematica/mathematica.hpp"
#include "part.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {

template<typename Frame>
Part<Frame>::Part(
    DegreesOfFreedom<Frame> const& degrees_of_freedom,
    Mass const& mass,
    Vector<Acceleration, Frame> const&
        gravitational_acceleration_to_be_applied_by_ksp)
    : degrees_of_freedom_(degrees_of_freedom),
      mass_(mass),
      gravitational_acceleration_to_be_applied_by_ksp_(
          gravitational_acceleration_to_be_applied_by_ksp) {}

template<typename Frame>
DegreesOfFreedom<Frame> const& Part<Frame>::degrees_of_freedom() const {
  return degrees_of_freedom_;
}

template<typename Frame>
Mass const& Part<Frame>::mass() const {
  return mass_;
}

template<typename Frame>
Vector<Acceleration, Frame> const&
    Part<Frame>::gravitational_acceleration_to_be_applied_by_ksp() const {
  return gravitational_acceleration_to_be_applied_by_ksp_;
}

template<typename Frame>
void Part<Frame>::WriteToMessage(
    not_null<serialization::Part*> const message) const {
  degrees_of_freedom_.WriteToMessage(message->mutable_degrees_of_freedom());
  mass_.WriteToMessage(message->mutable_mass());
  gravitational_acceleration_to_be_applied_by_ksp_.WriteToMessage(
      message->mutable_gravitational_acceleration_to_be_applied_by_ksp());
}

template<typename Frame>
Part<Frame> Part<Frame>::ReadFromMessage(serialization::Part const& message) {
  return Part(DegreesOfFreedom<Frame>::ReadFromMessage(
                  message.degrees_of_freedom()),
              Mass::ReadFromMessage(message.mass()),
              Vector<Acceleration, Frame>::ReadFromMessage(
                  message.gravitational_acceleration_to_be_applied_by_ksp()));
}

template<typename Frame>
std::ostream& operator<<(std::ostream& out, Part<Frame> const& part) {
  return out << "{"
      << part.degrees_of_freedom() << ", "
      << part.mass() << ", "
      << part.gravitational_acceleration_to_be_applied_by_ksp() << "}";
}

}  // namespace ksp_plugin
}  // namespace principia
