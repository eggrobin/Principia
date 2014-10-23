#pragma once

#include "geometry/affine_map.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "ksp_plugin/celestial.hpp"

using principia::geometry::AngularVelocity;

namespace principia {
namespace ksp_plugin {

struct Barycentre;

class RenderingFrame {
 public:
  // Returns a trajectory representing |actual_trajectory| as seen in the
  // rendering frame, realized in the |Barycentre| frame.
  virtual std::unique_ptr<Trajectory<Barycentre>> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const = 0;
 protected:
  // The rendering frame.  The position and rotation with respect to
  // |Barycentre| is implementation-defined and may depend on time.
  // The orientation should be the same as that of |Barycentre|.
  // NOTE(egg): "Apparent reference frame" is probably a Gallicism.
  struct Apparent;
  static Position<Apparent> const RenderingFrameOrigin() {
    return Position<Apparent>();
  }
};

class BodyCentredNonRotatingFrame : public RenderingFrame {
 public:
    
  explicit BodyCentredNonRotatingFrame(Celestial<Barycentre> const& body);

  std::unique_ptr<Trajectory<Barycentre>> const ApparentTrajectory(
      Trajectory<Barycentre> const& actual_trajectory) const override;

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
