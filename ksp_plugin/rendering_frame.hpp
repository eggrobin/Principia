#pragma once

#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "ksp_plugin/celestial.hpp"

using principia::geometry::AngularVelocity;

namespace principia {
namespace ksp_plugin {

// Functions for converting between a reference frame |MovingFrame| in motion
// with respect to |FixedFrame| and |FixedFrame|.
template<typename MovingFrame, typename FixedFrame>
class MovingFrameDefinition {
 public:
  virtual std::unique_ptr<Trajectory<MovingFrame>> const ApparentTrajectory(
      Trajectory<FixedFrame> const& actual_trajectory) const = 0;
  virtual Position<MovingFrame> const ToMovingFrame(Instant const& t) const = 0;
  virtual Position<MovingFrame> const ToMovingFrameAtCurrentTime(
      Position<FixedFrame> const& position) const = 0;
  virtual DegreesOfFreedom<MovingFrame> const ToMovingFrame(
      DegreesOfFreedom<FixedFrame> const& state,
      Instant const& t);
  virtual Rotation<FixedFrame, MovingFrame> RotationOfMovingFrame(
      Instant const& t) const = 0;
  virtual Rotation<FixedFrame, MovingFrame>
  RotationOfMovingFrameAtCurrentTime() const = 0;
};


struct Barycentre;

// The reference frame in which trajectories are rendered. Note that we may wish
// to render different trajectories in different referent frame at the same
// time, so that the type alone cannot specify the reference frame.
// Templates using this as a  template parametre should be wrapped together with
// a subclass of |RenderingFrameDefinition| in order not to be ambiguous.
struct RenderingFrame;
static Position<RenderingFrame> const kRenderingFrameOrigin;

// Reference frames used for rendering are always defined with respect to
// |Barycentre|.
using RenderingFrameDefinition = MovingFrameDefinition<RenderingFrame,
                                                       Barycentre> {};

class BodyCentredNonRotatingFrame : RenderingFrameDefinition {
 public:
  virtual std::unique_ptr<Trajectory<RenderingFrame>> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const = 0;
  virtual Position<RenderingFrame> const ToMovingFrame(
      Instant const& t) const = 0;
  virtual Position<RenderingFrame> const ToMovingFrameAtCurrentTime(
      Position<Barycentre> const& position) const = 0;
  virtual DegreesOfFreedom<RenderingFrame> const ToMovingFrame(
      DegreesOfFreedom<Barycentre> const& state,
      Instant const& t);
  virtual Rotation<Barycentre, RenderingFrame> RotationOfMovingFrame(
      Instant const& t) const = 0;
  virtual Rotation<Barycentre, RenderingFrame>
  RotationOfRenderingFrameAtCurrentTime() const = 0;

 private:
  Celestial<Barycentre> const& body_;
};

class BodyCentredRotatingWithSurface : RenderingFrame {
 public:
  BodyCentredRotatingWithSurface(Celestial<Barycentre> const& body,
                                 AngularVelocity<Barycentre> const& rotation_);

  std::unique_ptr<Trajectory<Barycentre>> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const override;

 private:
  Celestial<Barycentre> const& body_;
  AngularVelocity<Barycentre> const rotation_;
};

class BarycentricRotating : RenderingFrame {
 public:
  BarycentricRotating(Celestial<Barycentre> const& primary,
                      Celestial<Barycentre> const& secondary_);

  std::unique_ptr<Trajectory<Barycentre>> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const override;

 private:
  Celestial<Barycentre> const& primary_;
  Celestial<Barycentre> const& secondary_;
};

}  // namespace ksp_plugin
}  // namespace principia
