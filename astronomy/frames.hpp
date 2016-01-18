﻿
#pragma once

#include <functional>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using geometry::Frame;
using geometry::Instant;
using geometry::Position;
using quantities::si::ArcMinute;
using quantities::si::ArcSecond;
using quantities::si::Degree;

namespace astronomy {

// A reference frame with a basis.
// The frame is the International Celestial Reference Frame.
// The basis is defined from the orbit of the Earth at J2000.0 as follows:
// The xy plane is the plane of the Earth's orbit at J2000.0.
// The x axis is out along the ascending node of the instantaneous plane of the
// Earth's orbit and the Earth's mean equator at J2000.0.
// The z axis is perpendicular to the xy-plane in the directional (+ or -) sense
// of Earth's north pole at J2000.0.
// The basis is right-handed and orthonormal.
using ICRFJ2000Ecliptic =
    Frame<serialization::Frame::SolarSystemTag,
          serialization::Frame::ICRF_J2000_ECLIPTIC, true>;

// The xy plane is the plane of the Earth's mean equator at J2000.0.
// The x axis is out along the ascending node of the instantaneous plane of the
// Earth's orbit and the Earth's mean equator at J2000.0.
// The z axis is along the Earth's mean north pole at J2000.0.
// The basis is right-handed and orthonormal.
// Note that |ICRFJ2000Equator| and |ICRFJ2000Ecliptic| share their x axis.
using ICRFJ2000Equator =
    Frame<serialization::Frame::SolarSystemTag,
          serialization::Frame::ICRF_J2000_EQUATOR, true>;

// Rotation around the common x axis mapping equatorial coordinates to ecliptic
// coordinates.  The angle is the one defined by the XVIth General Assembly of
// the International Astronomical Union.
geometry::Rotation<ICRFJ2000Equator, ICRFJ2000Ecliptic> const
    kEquatorialToEcliptic =
        geometry::Rotation<ICRFJ2000Equator, ICRFJ2000Ecliptic>(
            23 * Degree + 26 * ArcMinute + 21.448 * ArcSecond,
            geometry::Bivector<double, ICRFJ2000Equator>({-1, 0, 0}));

geometry::Position<ICRFJ2000Ecliptic> const kSolarSystemBarycentreEcliptic;
geometry::Position<ICRFJ2000Equator> const kSolarSystemBarycentreEquator;

}  // namespace astronomy
}  // namespace principia
