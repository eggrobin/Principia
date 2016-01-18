
#pragma once

#include "base/not_null.hpp"
#include "journal/player.hpp"
#include "ksp_plugin/burn.hpp"
#include "ksp_plugin/interface.hpp"
#include "serialization/journal.pb.h"

struct Burn;
struct KSPPart;
struct NavigationFrameParameters;
struct NavigationManoeuvre;
struct QP;
struct WXYZ;
struct XYZ;
struct XYZSegment;

namespace principia {

namespace base {
template <typename Pointer> class not_null;
}  // namespace base
namespace interface {
struct LineAndIterator;
}  // namespace interface

using base::not_null;
using interface::Burn;
using interface::KSPPart;
using interface::LineAndIterator;
using interface::NavigationFrameParameters;
using interface::NavigationManoeuvre;
using interface::QP;
using interface::WXYZ;
using interface::XYZ;
using interface::XYZSegment;

namespace journal {

#include "journal/profiles.generated.h"

}  // namespace journal
}  // namespace principia
