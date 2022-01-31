module;

#include <utility>

#include "physics/continuous_trajectory.hpp"
#include "physics/rotating_body.hpp"

module principia.ksp_plugin.celestial;

namespace principia {

using base::not_null;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
using physics::RotatingBody;

namespace ksp_plugin {

Celestial::Celestial(not_null<RotatingBody<Barycentric> const*> body)
    : body_(std::move(body)) {}

bool Celestial::is_initialized() const {
  return trajectory_ != nullptr;
}

void Celestial::set_trajectory(
    not_null<ContinuousTrajectory<Barycentric> const*> const trajectory) {
  CHECK(!is_initialized());
  trajectory_ = trajectory;
}

ContinuousTrajectory<Barycentric> const& Celestial::trajectory() const {
  CHECK(is_initialized());
  return *trajectory_;
}

DegreesOfFreedom<Barycentric> Celestial::current_degrees_of_freedom(
    Instant const& current_time) const {
  CHECK(is_initialized());
  return trajectory().EvaluateDegreesOfFreedom(current_time);
}

Position<Barycentric> Celestial::current_position(
    Instant const& current_time) const {
  CHECK(is_initialized());
  return trajectory().EvaluatePosition(current_time);
}

Velocity<Barycentric> Celestial::current_velocity(
    Instant const& current_time) const {
  CHECK(is_initialized());
  return trajectory().EvaluateVelocity(current_time);
}

not_null<RotatingBody<Barycentric> const*> Celestial::body() const {
  return body_;
}

bool Celestial::has_parent() const {
  return parent_ != nullptr;
}

Celestial const* Celestial::parent() const {
  return parent_;
}

void Celestial::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

}  // namespace ksp_plugin
}  // namespace principia
