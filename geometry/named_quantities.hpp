﻿
#pragma once

#include <type_traits>

#include "geometry/affine_map.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/point.hpp"
#include "geometry/symmetric_bilinear_form.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {

using Instant = Point<quantities::Time>;
template<typename Frame>
using Displacement = Vector<quantities::Length, Frame>;
template<typename Frame>
using Position = Point<Displacement<Frame>>;
template<typename Frame>
using Velocity = Vector<quantities::Speed, Frame>;

template<typename Frame>
using AngularVelocity = Bivector<quantities::AngularFrequency, Frame>;

// An arbitrary rigid transformation.  Simultaneous positions between two frames
// are always related by such a transformation.
template<typename FromFrame, typename ToFrame>
using RigidTransformation =
    AffineMap<FromFrame, ToFrame, quantities::Length, OrthogonalMap>;

template<typename Frame>
using InertiaTensor =
    SymmetricBilinearForm<quantities::MomentOfInertia, Frame, Bivector>;

// IEEE 754:2008 nextUp and nextDown for Instants.
// We would like to avoid the terms “up” and “down” when referring to the
// passage of time.  We avoid the term “next” in one direction because of the
// confusability with |std::nextafter|, which has different semantics, and in
// the other because of the awkwardness of the phrase “next before”.
// Defined inline for want of a way to alias functions in C++.
constexpr Instant JustAfter(Instant const t) { return NextUp(t); }
constexpr Instant JustBefore(Instant const t) { return NextDown(t); }

namespace internal_point {
// We must declare this in the internal namespace where Point is defined so that
// it is found by ADL.
std::ostream& operator<<(std::ostream& os, const Instant& t);
}  // namespace internal_point
}  // namespace geometry
}  // namespace principia
