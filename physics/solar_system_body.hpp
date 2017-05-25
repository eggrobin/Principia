﻿
#pragma once

#include "physics/solar_system.hpp"

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/time_scales.hpp"
#include "base/map_util.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/hierarchical_system.hpp"
#include "physics/massive_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/parser.hpp"
#include "quantities/si.hpp"
#include "serialization/astronomy.pb.h"

namespace principia {
namespace physics {
namespace internal_solar_system {

using astronomy::J2000;
using astronomy::ParseTT;
using base::Contains;
using base::dynamic_cast_not_null;
using base::FindOrDie;
using base::make_not_null_unique;
using geometry::Bivector;
using geometry::Frame;
using geometry::Instant;
using geometry::Position;
using geometry::RadiusLatitudeLongitude;
using geometry::Vector;
using geometry::Velocity;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::DebugString;
using quantities::Length;
using quantities::ParseQuantity;
using quantities::Speed;
using quantities::si::Radian;
using quantities::si::Second;

inline serialization::GravityModel ParseGravityModel(
    std::experimental::filesystem::path const& gravity_model_filename) {
  serialization::SolarSystemFile gravity_model;
  std::ifstream gravity_model_ifstream(gravity_model_filename);
  CHECK(gravity_model_ifstream.good());
  google::protobuf::io::IstreamInputStream gravity_model_zcs(
                                               &gravity_model_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&gravity_model_zcs,
                                            &gravity_model));
  CHECK(gravity_model.has_gravity_model());
  return gravity_model.gravity_model();
}

inline serialization::InitialState ParseInitialState(
    std::experimental::filesystem::path const& initial_state_filename) {
  serialization::SolarSystemFile initial_state;
  std::ifstream initial_state_ifstream(initial_state_filename);
  CHECK(initial_state_ifstream.good());
  google::protobuf::io::IstreamInputStream initial_state_zcs(
                                               &initial_state_ifstream);
  CHECK(google::protobuf::TextFormat::Parse(&initial_state_zcs,
                                            &initial_state));
  CHECK(initial_state.has_initial_state());
  return initial_state.initial_state();
}

template<typename Frame>
SolarSystem<Frame>::SolarSystem(
    std::experimental::filesystem::path const& gravity_model_filename,
    std::experimental::filesystem::path const& initial_state_filename)
    : SolarSystem(ParseGravityModel(gravity_model_filename),
                  ParseInitialState(initial_state_filename)) {}

template<typename Frame>
SolarSystem<Frame>::SolarSystem(
    serialization::GravityModel const& gravity_model,
    serialization::InitialState const& initial_state)
    : gravity_model_(gravity_model),
      initial_state_(initial_state) {
  gravity_model_.CheckInitialized();
  initial_state_.CheckInitialized();

  CheckFrame(gravity_model_);
  CheckFrame(initial_state_);

  // Store the data in maps keyed by body name.
  for (auto& body : *gravity_model_.mutable_body()) {
    bool inserted;
    std::tie(std::ignore, inserted) =
        gravity_model_map_.insert(std::make_pair(body.name(), &body));
    CHECK(inserted) << body.name();
  }
  if (initial_state_.has_cartesian()) {
    for (auto const& body : initial_state_.cartesian().body()) {
      bool inserted;
      std::tie(std::ignore, inserted) =
          cartesian_initial_state_map_.emplace(body.name(), &body);
      CHECK(inserted) << body.name();
    }

    // Check that the maps are consistent.
    auto it1 = gravity_model_map_.begin();
    auto it2 = cartesian_initial_state_map_.begin();
    for (; it1 != gravity_model_map_.end() &&
           it2 != cartesian_initial_state_map_.end();
         ++it1, ++it2) {
      CHECK_EQ(it1->first, it2->first);
      names_.push_back(it1->first);
    }
    CHECK(it1 == gravity_model_map_.end()) << it1->first;
    CHECK(it2 == cartesian_initial_state_map_.end()) << it2->first;
  } else {
    for (auto& body : *initial_state_.mutable_keplerian()->mutable_body()) {
      bool inserted;
      std::tie(std::ignore, inserted) =
          keplerian_initial_state_map_.emplace(body.name(), &body);
      CHECK(inserted) << body.name();
    }

    // Check that the maps are consistent.
    auto it1 = gravity_model_map_.begin();
    auto it2 = keplerian_initial_state_map_.begin();
    for (; it1 != gravity_model_map_.end() &&
           it2 != keplerian_initial_state_map_.end();
         ++it1, ++it2) {
      CHECK_EQ(it1->first, it2->first);
      names_.push_back(it1->first);
    }
    CHECK(it1 == gravity_model_map_.end()) << it1->first;
    CHECK(it2 == keplerian_initial_state_map_.end()) << it2->first;
  }

  epoch_ = ParseTT(initial_state_.epoch());

  // Call these two functions to parse all the data, so that errors are detected
  // at initialization.  Drop their results on the floor.
  MakeAllMassiveBodies();
  MakeAllDegreesOfFreedom();
}

template<typename Frame>
std::unique_ptr<Ephemeris<Frame>> SolarSystem<Frame>::MakeEphemeris(
    Length const& fitting_tolerance,
    typename Ephemeris<Frame>::FixedStepParameters const& parameters) {
  return std::make_unique<Ephemeris<Frame>>(MakeAllMassiveBodies(),
                                            MakeAllDegreesOfFreedom(),
                                            epoch_,
                                            fitting_tolerance,
                                            parameters);
}

template<typename Frame>
Instant const& SolarSystem<Frame>::epoch() const {
  return epoch_;
}

template<typename Frame>
std::vector<std::string> const& SolarSystem<Frame>::names() const {
  return names_;
}

template<typename Frame>
int SolarSystem<Frame>::index(std::string const& name) const {
  auto const it = std::equal_range(names_.begin(), names_.end(), name);
  return it.first - names_.begin();
}


template<typename Frame>
DegreesOfFreedom<Frame> SolarSystem<Frame>::degrees_of_freedom(
    std::string const& name) const {
  return MakeDegreesOfFreedom(*cartesian_initial_state_map_.at(name));
}

template<typename Frame>
GravitationalParameter SolarSystem<Frame>::gravitational_parameter(
    std::string const& name) const {
  return MakeMassiveBody(*gravity_model_map_.at(name))->
             gravitational_parameter();
}

template<typename Frame>
Length SolarSystem<Frame>::mean_radius(std::string const& name) const {
  return MakeMassiveBody(*gravity_model_map_.at(name))->mean_radius();
}

template<typename Frame>
not_null<MassiveBody const*> SolarSystem<Frame>::massive_body(
    Ephemeris<Frame> const & ephemeris,
    std::string const & name) const {
  return ephemeris.bodies()[index(name)];
}

template<typename Frame>
not_null<RotatingBody<Frame> const*> SolarSystem<Frame>::rotating_body(
    Ephemeris<Frame> const& ephemeris,
    std::string const& name) const {
  CHECK(gravity_model_message(name).has_mean_radius());
  return dynamic_cast_not_null<RotatingBody<Frame> const*>(
      massive_body(ephemeris, name));
}

template<typename Frame>
ContinuousTrajectory<Frame> const& SolarSystem<Frame>::trajectory(
    Ephemeris<Frame> const & ephemeris,
    std::string const & name) const {
  MassiveBody const* const body = ephemeris.bodies()[index(name)];
  return *ephemeris.trajectory(body);
}

template<typename Frame>
serialization::GravityModel::Body const&
SolarSystem<Frame>::gravity_model_message(std::string const& name) const {
  return *FindOrDie(gravity_model_map_, name);
}

template<typename Frame>
serialization::InitialState::Cartesian::Body const&
SolarSystem<Frame>::cartesian_initial_state_message(
    std::string const& name) const {
  return *FindOrDie(cartesian_initial_state_map_, name);
}

template<typename Frame>
serialization::InitialState::Keplerian::Body const&
SolarSystem<Frame>::keplerian_initial_state_message(
    std::string const& name) const {
  return *FindOrDie(keplerian_initial_state_map_, name);
}

template<typename Frame>
DegreesOfFreedom<Frame> SolarSystem<Frame>::MakeDegreesOfFreedom(
    serialization::InitialState::Cartesian::Body const& body) {
  Position<Frame> const
      position = Frame::origin +
                 Vector<Length, Frame>({ParseQuantity<Length>(body.x()),
                                        ParseQuantity<Length>(body.y()),
                                        ParseQuantity<Length>(body.z())});
  Velocity<Frame> const
      velocity(Vector<Speed, Frame>({ParseQuantity<Speed>(body.vx()),
                                     ParseQuantity<Speed>(body.vy()),
                                     ParseQuantity<Speed>(body.vz())}));
  return DegreesOfFreedom<Frame>(position, velocity);
}

template<typename Frame>
KeplerianElements<Frame> SolarSystem<Frame>::MakeKeplerianElements(
    serialization::InitialState::Keplerian::Body::Elements const& elements) {
  KeplerianElements<Frame> result;
  result.eccentricity = elements.eccentricity();
  CHECK_NE(elements.has_semimajor_axis(), elements.has_mean_motion());
  if (elements.has_semimajor_axis()) {
    result.semimajor_axis = ParseQuantity<Length>(elements.semimajor_axis());
  }
  if (elements.has_mean_motion()) {
    result.mean_motion =
        ParseQuantity<AngularFrequency>(elements.mean_motion());
  }
  result.inclination = ParseQuantity<Angle>(elements.inclination());
  result.longitude_of_ascending_node =
      ParseQuantity<Angle>(elements.longitude_of_ascending_node());
  result.argument_of_periapsis =
      ParseQuantity<Angle>(elements.argument_of_periapsis());
  result.mean_anomaly = ParseQuantity<Angle>(elements.mean_anomaly());
  return result;
}

template<typename Frame>
not_null<std::unique_ptr<MassiveBody>> SolarSystem<Frame>::MakeMassiveBody(
    serialization::GravityModel::Body const& body) {
  Check(body);
  auto const massive_body_parameters = MakeMassiveBodyParameters(body);
  if (body.has_mean_radius()) {
    auto const rotating_body_parameters = MakeRotatingBodyParameters(body);
    if (body.has_j2()) {
      auto const oblate_body_parameters = MakeOblateBodyParameters(body);
      return std::make_unique<OblateBody<Frame>>(*massive_body_parameters,
                                                 *rotating_body_parameters,
                                                 *oblate_body_parameters);
    } else {
      return std::make_unique<RotatingBody<Frame>>(*massive_body_parameters,
                                                   *rotating_body_parameters);
    }
  } else {
    return std::make_unique<MassiveBody>(*massive_body_parameters);
  }
}

template<typename Frame>
not_null<std::unique_ptr<RotatingBody<Frame>>>
SolarSystem<Frame>::MakeRotatingBody(
    serialization::GravityModel::Body const& body) {
  Check(body);
  CHECK(body.has_mean_radius());
  auto const massive_body_parameters = MakeMassiveBodyParameters(body);
  auto const rotating_body_parameters = MakeRotatingBodyParameters(body);
  if (body.has_j2()) {
    auto const oblate_body_parameters = MakeOblateBodyParameters(body);
    return std::make_unique<OblateBody<Frame>>(*massive_body_parameters,
                                               *rotating_body_parameters,
                                               *oblate_body_parameters);
  } else {
    return std::make_unique<RotatingBody<Frame>>(*massive_body_parameters,
                                                 *rotating_body_parameters);
  }
}

template<typename Frame>
not_null<std::unique_ptr<OblateBody<Frame>>>
SolarSystem<Frame>::MakeOblateBody(
    serialization::GravityModel::Body const& body) {
  Check(body);
  CHECK(body.has_mean_radius());
  CHECK(body.has_j2());
  auto const massive_body_parameters = MakeMassiveBodyParameters(body);
  auto const rotating_body_parameters = MakeRotatingBodyParameters(body);
  auto const oblate_body_parameters = MakeOblateBodyParameters(body);
  return std::make_unique<OblateBody<Frame>>(*massive_body_parameters,
                                             *rotating_body_parameters,
                                             *oblate_body_parameters);
}

template<typename Frame>
not_null<std::unique_ptr<HierarchicalSystem<Frame>>>
SolarSystem<Frame>::MakeHierarchicalSystem() {
  // First, construct all the bodies and find the primary body of the system.
  std::string primary;
  std::map<std::string,
            not_null<std::unique_ptr<MassiveBody const>>> owned_bodies;
  std::map<std::string, not_null<MassiveBody const*>> unowned_bodies;
  for (auto const& pair : keplerian_initial_state_map_) {
    const auto& name = pair.first;
    serialization::InitialState::Keplerian::Body* const body = pair.second;
    CHECK_EQ(body->has_parent(), body->has_elements()) << name;
    if (!body->has_parent()) {
      CHECK(primary.empty()) << name;
      primary = name;
    }
    auto owned_body = MakeMassiveBody(gravity_model_message(name));
    unowned_bodies.emplace(name, owned_body.get());
    owned_bodies.emplace(name, std::move(owned_body));
  }

  // Construct a hierarchical system rooted at the primary and add the other
  // bodies layer by layer.
  auto hierarchical_system = make_not_null_unique<HierarchicalSystem<Frame>>(
      std::move(FindOrDie(owned_bodies, primary)));
  std::set<std::string> previous_layer = {primary};
  std::set<std::string> current_layer;
  do {
    for (auto const& pair : keplerian_initial_state_map_) {
      const auto& name = pair.first;
      serialization::InitialState::Keplerian::Body* const body = pair.second;
      if (Contains(previous_layer, body->parent())) {
        current_layer.insert(name);
        KeplerianElements<Frame> const elements =
            MakeKeplerianElements(body->elements());
        hierarchical_system->Add(std::move(FindOrDie(owned_bodies, name)),
                                 FindOrDie(unowned_bodies, body->parent()),
                                 elements);
      }
    }
    previous_layer = current_layer;
    current_layer.clear();
  } while (!previous_layer.empty());

  return hierarchical_system;
}

template<typename Frame>
void SolarSystem<Frame>::RemoveMassiveBody(std::string const& name) {
  for (int i = 0; i < names_.size(); ++i) {
    if (names_[i] == name) {
      names_.erase(names_.begin() + i);
      cartesian_initial_state_map_.erase(name);
      gravity_model_map_.erase(name);
      return;
    }
  }
  LOG(FATAL) << name << " does not exist";
}

template<typename Frame>
void SolarSystem<Frame>::RemoveOblateness(std::string const& name) {
  auto const it = gravity_model_map_.find(name);
  CHECK(it != gravity_model_map_.end()) << name << " does not exist";
  serialization::GravityModel::Body* body = it->second;
  body->clear_j2();
  body->clear_reference_radius();
}

template<typename Frame>
void SolarSystem<Frame>::ReplaceElements(
    std::string const& name,
    KeplerianElements<Frame> const& elements) {
  auto* const body_elements =
      FindOrDie(keplerian_initial_state_map_, name)->mutable_elements();
  body_elements->set_eccentricity(*elements.eccentricity);
  if (elements.semimajor_axis) {
    body_elements->set_semimajor_axis(DebugString(*elements.semimajor_axis));
  }
  if (elements.mean_motion) {
    body_elements->set_mean_motion(DebugString(*elements.mean_motion));
  }
  body_elements->set_inclination(DebugString(elements.inclination));
  body_elements->set_longitude_of_ascending_node(
      DebugString(elements.longitude_of_ascending_node));
  body_elements->set_argument_of_periapsis(
      DebugString(*elements.argument_of_periapsis));
  body_elements->set_mean_anomaly(DebugString(*elements.mean_anomaly));
}

template<typename Frame>
void SolarSystem<Frame>::Check(serialization::GravityModel::Body const& body) {
  CHECK(body.has_name());
  CHECK(body.has_gravitational_parameter()) << body.name();
  CHECK_EQ(body.has_reference_instant(), body.has_mean_radius()) << body.name();
  CHECK_EQ(body.has_reference_instant(),
           body.has_axis_right_ascension()) << body.name();
  CHECK_EQ(body.has_reference_instant(),
           body.has_axis_declination()) << body.name();
  CHECK_EQ(body.has_reference_instant(),
           body.has_reference_angle()) << body.name();
  CHECK_EQ(body.has_reference_instant(),
           body.has_angular_frequency()) << body.name();
  CHECK_EQ(body.has_j2(), body.has_reference_radius()) << body.name();
}

template<typename Frame>
not_null<std::unique_ptr<MassiveBody::Parameters>>
SolarSystem<Frame>::MakeMassiveBodyParameters(
    serialization::GravityModel::Body const& body) {
  return make_not_null_unique<MassiveBody::Parameters>(
      body.name(),
      ParseQuantity<GravitationalParameter>(body.gravitational_parameter()));
}

template<typename Frame>
not_null<std::unique_ptr<typename RotatingBody<Frame>::Parameters>>
SolarSystem<Frame>::MakeRotatingBodyParameters(
    serialization::GravityModel::Body const& body) {
  return make_not_null_unique<typename RotatingBody<Frame>::Parameters>(
      ParseQuantity<Length>(body.mean_radius()),
      ParseQuantity<Angle>(body.reference_angle()),
      ParseTT(body.reference_instant()),
      ParseQuantity<AngularFrequency>(body.angular_frequency()),
      ParseQuantity<Angle>(body.axis_right_ascension()),
      ParseQuantity<Angle>(body.axis_declination()));
}

template<typename Frame>
not_null<std::unique_ptr<typename OblateBody<Frame>::Parameters>>
SolarSystem<Frame>::MakeOblateBodyParameters(
    serialization::GravityModel::Body const& body) {
  return make_not_null_unique<typename OblateBody<Frame>::Parameters>(
      body.j2(),
      ParseQuantity<Length>(body.reference_radius()));
}

template<typename Frame>
std::vector<not_null<std::unique_ptr<MassiveBody const>>>
SolarSystem<Frame>::MakeAllMassiveBodies() {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  for (auto const& pair : gravity_model_map_) {
    serialization::GravityModel::Body const* const body = pair.second;
    bodies.emplace_back(MakeMassiveBody(*body));
  }
  return bodies;
}

template<typename Frame>
std::vector<DegreesOfFreedom<Frame>>
SolarSystem<Frame>::MakeAllDegreesOfFreedom() {
  std::vector<DegreesOfFreedom<Frame>> degrees_of_freedom;
  if (!cartesian_initial_state_map_.empty()) {
    for (auto const& pair : cartesian_initial_state_map_) {
      serialization::InitialState::Cartesian::Body const* const body =
          pair.second;
      degrees_of_freedom.push_back(MakeDegreesOfFreedom(*body));
    }
  }
  if (!keplerian_initial_state_map_.empty()) {
    auto const hierarchical_system = MakeHierarchicalSystem();

    // Construct a barycentric system and fill a map from body name to degrees
    // of freedom.
    typename HierarchicalSystem<Frame>::BarycentricSystem const
        barycentric_system = hierarchical_system->ConsumeBarycentricSystem();
    std::map<std::string, DegreesOfFreedom<Frame>> name_to_degrees_of_freedom;
    for (int i = 0; i < barycentric_system.bodies.size(); ++i) {
      auto const& body = barycentric_system.bodies[i];
      auto const& degrees_of_freedom = barycentric_system.degrees_of_freedom[i];
      name_to_degrees_of_freedom.emplace(body->name(), degrees_of_freedom);
    }

    // Obtain the final result with degrees of freedom in name order.
    for (auto const& pair : name_to_degrees_of_freedom) {
      degrees_of_freedom.push_back(pair.second);
    }
  }
  return degrees_of_freedom;
}

template<typename Frame>
template<typename Message>
void SolarSystem<Frame>::CheckFrame(Message const& message) {
  if (message.has_solar_system_frame()) {
    CHECK_EQ(Frame::tag, message.solar_system_frame());
  }
  if (message.has_plugin_frame()) {
    CHECK_EQ(Frame::tag, message.plugin_frame());
  }
}

}  // namespace internal_solar_system
}  // namespace physics
}  // namespace principia
