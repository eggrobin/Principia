#pragma once

#include "ksp_plugin/manœuvre.hpp"

#include <cmath>

#include "quantities/elementary_functions.hpp"

namespace principia {

using quantities::Sqrt;

namespace ksp_plugin {

inline Manœuvre::Manœuvre(Force thrust,
                          Mass initial_mass,
                          Speed effective_exhaust_velocity,
                          Vector<double, Frenet> direction)
    : thrust_(thrust),
      initial_mass_(initial_mass),
      effective_exhaust_velocity_(effective_exhaust_velocity),
      frenet_initial_direction_(direction) {}

inline Instant Manœuvre::initial_time() const {
  CHECK(initial_time_ != nullptr);
  return *initial_time_;
}

inline Instant Manœuvre::time_of_half_Δv() const {
  return initial_time() + time_to_half_Δv();
}

inline Instant Manœuvre::final_time() const {
  return initial_time() + duration();
}

inline void Manœuvre::set_initial_time(Instant const & initial_time) {
  initial_time_ = std::make_unique<Instant>(initial_time);
}

inline void Manœuvre::set_time_of_half_Δv(Instant const& time_of_half_Δv) {
  set_initial_time(time_of_half_Δv - time_to_half_Δv());
}

inline Time Manœuvre::duration() const {
  CHECK(duration_ != nullptr);
  return *duration_;
}

inline void Manœuvre::set_duration(Time const & duration) {
  duration_ = std::make_unique<Time>(duration);
}

inline void Manœuvre::set_Δv(Speed const& Δv) {
  set_duration(initial_mass() * effective_exhaust_velocity() *
               (1 - std::exp(-Δv / effective_exhaust_velocity())) / thrust());
}

inline Speed Manœuvre::Δv() const {
  // Циолко́вский's equation.
  return effective_exhaust_velocity() * std::log(initial_mass() / final_mass());
}

inline Vector<double, Frenet> Manœuvre::direction() const {
  return frenet_initial_direction_;
}

inline void Manœuvre::SetFrenetFrame(
    Rotation<Frenet, Barycentric> const& rotation) {
  inertial_direction_ = rotation(frenet_initial_direction_);
}

inline Speed Manœuvre::effective_exhaust_velocity() const {
  return effective_exhaust_velocity_;
}

inline Force Manœuvre::thrust() const {
  return thrust_;
}

inline Mass Manœuvre::initial_mass() const {
  return initial_mass_;
}

inline Variation<Mass> Manœuvre::mass_flow() const {
  return thrust() / effective_exhaust_velocity();
}

inline Mass Manœuvre::final_mass() const {
  return initial_mass() - mass_flow() * duration();
}

inline Time Manœuvre::time_to_half_Δv() const {
  return effective_exhaust_velocity() * initial_mass() *
         (1 - std::sqrt(final_mass() / initial_mass())) / thrust();
}

inline DiscreteTrajectory<Barycentric>::IntrinsicAcceleration
Manœuvre::acceleration() const {
  return [
    direction = *inertial_direction_,
    initial_time = this->initial_time(),
    final_time = this->final_time(),
    thrust = this->thrust(),
    initial_mass = this->initial_mass(),
    mass_flow = this->mass_flow()
  ](Instant const& time) -> Vector<Acceleration, Barycentric> {
    if (time >= initial_time && time <= final_time) {
      return direction * thrust /
             (initial_mass - (time - initial_time) * mass_flow);
    } else {
      return Vector<Acceleration, Barycentric>();
    }
  };
}

}  // namespace ksp_plugin
}  // namespace principia
