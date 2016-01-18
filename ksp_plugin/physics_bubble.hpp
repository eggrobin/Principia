﻿#pragma once

#include <stddef.h>
#include <experimental/optional>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/rotation.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/forkable.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "serialization/ksp_plugin.pb.h"

namespace principia {
namespace ksp_plugin {
class Vessel;
}  // namespace ksp_plugin
}  // namespace principia

namespace principia {

namespace geometry {
template <typename Scalar, typename Frame, int rank> class Multivector;
template <typename Vector> class Point;
}  // namespace geometry
namespace ksp_plugin {
class Celestial;
}  // namespace ksp_plugin
namespace serialization {
class PhysicsBubble;
}  // namespace serialization

using base::not_null;
using geometry::OrthogonalMap;
using physics::DegreesOfFreedom;
using physics::RelativeDegreesOfFreedom;

namespace ksp_plugin {

class PhysicsBubble {
 public:
  using BarycentricToWorldSun = geometry::OrthogonalMap<Barycentric, WorldSun>;

  PhysicsBubble();
  ~PhysicsBubble() = default;

  // Creates |next_| if it is null.  Adds the |vessel| to |next_->vessels| with
  // a list of pointers to the Parts in |parts|.  Merges |parts| into
  // |next_->parts|.  The |vessel| must not already be in |next_->vessels|.
  // |parts| must not contain a |PartId| already in |next_->parts|.
  void AddVesselToNext(not_null<Vessel*> const vessel,
                       std::vector<IdAndOwnedPart> parts);

  // If |next_| is not null, computes the world centre of mass, trajectory
  // (including intrinsic acceleration) of |*next_|. Moves |next_| into
  // |current_|.  The trajectory of the centre of mass is reset to a single
  // point at |current_time| if the composition of the bubble changes.
  // TODO(phl): Document the parameters!
  void Prepare(BarycentricToWorldSun const& barycentric_to_world_sun,
               Instant const& current_time,
               Instant const& next_time);

  // Computes and returns |current_->displacement_correction|.  This is the
  // |World| shift to be applied to the bubble in order for it to be in the
  // correct position.
  Displacement<World> DisplacementCorrection(
      BarycentricToWorldSun const& barycentric_to_world_sun,
      Celestial const& reference_celestial,
      Position<World> const& reference_celestial_world_position) const;

  // Computes and returns |current_->velocity_correction|.  This is the |World|
  // shift to be applied to the physics bubble in order for it to have the
  // correct velocity.
  Velocity<World> VelocityCorrection(
      BarycentricToWorldSun const& barycentric_to_world_sun,
      Celestial const& reference_celestial) const;

  // Returns |current_ == nullptr|.
  bool empty() const;

  // Returns 0 if |empty()|, 1 otherwise.
  std::size_t count() const;

  // Returns |current_->vessels.size()|, or 0 if |empty()|.
  std::size_t number_of_vessels() const;

  // Returns true if, and only if, |vessel| is in |current_->vessels|.
  // |current_| may be null, in that case, returns false.
  bool contains(not_null<Vessel*> const vessel) const;

  // Selectors for the data in |current_|.
  std::vector<not_null<Vessel*>> vessels() const;
  RelativeDegreesOfFreedom<Barycentric> const& from_centre_of_mass(
      not_null<Vessel const*> const vessel) const;
  Ephemeris<Barycentric>::IntrinsicAcceleration const&
  centre_of_mass_intrinsic_acceleration() const;
  DiscreteTrajectory<Barycentric> const& centre_of_mass_trajectory() const;
  not_null<DiscreteTrajectory<Barycentric>*>
  mutable_centre_of_mass_trajectory() const;

  void WriteToMessage(
      std::function<std::string(not_null<Vessel const*>)> const guid,
      not_null<serialization::PhysicsBubble*> const message) const;
  static not_null<std::unique_ptr<PhysicsBubble>> ReadFromMessage(
      std::function<not_null<Vessel*>(std::string)> const vessel,
      serialization::PhysicsBubble const& message);

 private:
  using PartCorrespondence = std::pair<not_null<Part<World>*>,
                                       not_null<Part<World>*>>;

  struct PreliminaryState {
    PreliminaryState();
    std::map<not_null<Vessel*> const,
             // NOTE(Norgg) TODO(Egg) Removed const from vector,
             // custom allocator?
             std::vector<not_null<Part<World>*>>> vessels;
    PartIdToOwnedPart parts;
  };

  struct FullState : public PreliminaryState {
    explicit FullState(PreliminaryState preliminary_state);

    std::experimental::optional<DegreesOfFreedom<World>> centre_of_mass;
    Ephemeris<Barycentric>::IntrinsicAcceleration
        centre_of_mass_intrinsic_acceleration;
    std::unique_ptr<
        DiscreteTrajectory<Barycentric>> centre_of_mass_trajectory;
    std::experimental::optional<
        std::map<not_null<Vessel const*> const,
                 RelativeDegreesOfFreedom<Barycentric>>> from_centre_of_mass;
    std::experimental::optional<Displacement<World>> displacement_correction;
    std::experimental::optional<Velocity<World>> velocity_correction;
  };

  // Computes the world degrees of freedom of the centre of mass of
  // |next| using the contents of |next->parts|.
  void ComputeNextCentreOfMassWorldDegreesOfFreedom(
      not_null<FullState*> const next);

  // Computes |next->displacements_from_centre_of_mass| and
  // |next->velocities_from_centre_of_mass|.
  void ComputeNextVesselOffsets(
      BarycentricToWorldSun const& barycentric_to_world_sun,
      not_null<FullState*> const next);

  // Creates |next->centre_of_mass_trajectory| and appends to it the barycentre
  // of the degrees of freedom of the vessels in |next->vessels|.  There is no
  // intrinsic acceleration.
  void RestartNext(Instant const& current_time,
                   not_null<FullState*> const next);

  // Returns the parts common to |current_| and |next|.  The returned vector
  // contains pair of pointers to parts (current_part, next_part) for all parts
  // common to the two bubbles
  std::vector<PhysicsBubble::PartCorrespondence> ComputeCommonParts(
      FullState const& next);

  // Returns the intrinsic acceleration measured on the parts that are common to
  // the current and next bubbles.
  Vector<Acceleration, World> IntrinsicAcceleration(
      Instant const& current_time,
      Instant const& next_time,
      std::vector<PartCorrespondence>const& common_parts);

  // Given the vector of common parts, constructs
  // |next->centre_of_mass_trajectory| and appends degrees of freedom at
  // |current_time| that conserve the degrees of freedom of the centre of mass
  // of the parts in |common_parts|.
  void Shift(BarycentricToWorldSun const& barycentric_to_world_sun,
             Instant const& current_time,
             std::vector<PartCorrespondence> const& common_parts,
             not_null<FullState*> const next);

  std::unique_ptr<FullState> current_;
  // The following member is only accessed by |AddVesselToNext| and at the
  // beginning of |Prepare|.
  std::unique_ptr<PreliminaryState> next_;

  MasslessBody const body_;
};

}  // namespace ksp_plugin
}  // namespace principia
