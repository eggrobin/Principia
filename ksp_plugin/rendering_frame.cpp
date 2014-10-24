#include "ksp_plugin/rendering_frame.hpp"

#include "geometry/affine_map.hpp"
#include "geometry/rotation.hpp"
#include "ksp_plugin/celestial.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::AffineMap;
using principia::geometry::Rotation;
using principia::quantities::Time;

namespace principia {
namespace ksp_plugin {

Position<RenderingFrame::Apparent> const RenderingFrame::kRenderingFrameOrigin =
    Position<Apparent>;

BodyCentredNonRotatingFrame::BodyCentredNonRotatingFrame(
    Celestial<Barycentre> const& body) : body_(body) {};

std::unique_ptr<Trajectory<Barycentre>> const
BodyCentredNonRotatingFrame::ApparentTrajectory(
    Trajectory<Barycentre> const& actual_trajectory) const {
  std::unique_ptr<Trajectory<Barycentre>> result =
      std::make_unique<Trajectory<Barycentre>>(actual_trajectory.body());
  Trajectory<Barycentre>::Timeline const& actual_timeline =
      actual_trajectory.timeline();
  Trajectory<Barycentre>::Timeline const& reference_body_timeline =
      body_.history().timeline();
  Trajectory<Barycentre>::Timeline::const_iterator it_in_reference =
      reference_body_timeline.lower_bound(actual_timeline.begin()->first);
  DegreesOfFreedom<Barycentre> const& current_reference_state =
      body_.prolongation().timeline().rbegin()->second;
  CHECK(it_in_reference != reference_body_timeline.end());
  // TODO(egg): this is inelegant and probably a bit slow.  We need an
  // |IdentityMap| type.
  AffineMap<Barycentre, Apparent, Length, Rotation> const
  apparent_position_at_current_time =
      AffineMap<Apparent, Barycentre, Length, Rotation>(
          RenderingFrameOrigin(),
          current_reference_state.position,
          Rotation<Apparent, Barycentre>::Identity());
  for (auto const& pair : actual_timeline) {
    Instant const& t = pair.first;
    DegreesOfFreedom<Barycentre> const& actual_state = pair.second;
    while (it_in_reference->first < t) {
      ++it_in_reference;
    }
    if (it_in_reference->first == t) {
      DegreesOfFreedom<Barycentre> const& reference_state =
          it_in_reference->second;
      AffineMap<Barycentre, Apparent, Length, Rotation> apparent_position_at_t =
          AffineMap<Barycentre, Apparent, Length, Rotation>(
              reference_state.position,
              RenderingFrameOrigin(),
              Rotation<Barycentre, Apparent>::Identity());
      // TODO(egg): We should have a vector space structure on
      // |DegreesOfFreedom<Fries>|.
      result->Append(
          t,
          {apparent_position_at_current_time.Inverse()(
               apparent_position_at_t(actual_state.position)),
           actual_state.velocity - reference_state.velocity});
    }
  }
  return std::move(result);
}

}  // namespace ksp_plugin
}  // namespace principia
