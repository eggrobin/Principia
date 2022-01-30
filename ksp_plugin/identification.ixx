module;

#include <limits>
#include <map>
#include <set>
#include <string>

#include "base/not_null.hpp"

export module principia.ksp_plugin.identification;

namespace principia {

using base::not_null;

namespace ksp_plugin {

class Part;
class Vessel;

// The GUID of a vessel, obtained by |v.id.ToString()| in C#. We use this as a
// key in an |std::map|.
export using GUID = std::string;

// Corresponds to KSP's |Part.flightID|, *not* to |Part.uid|.  C#'s |uint|
// corresponds to |uint32_t|.
export using PartId = std::uint32_t;

// Comparator by PartId.  Useful for ensuring a consistent ordering in sets of
// pointers to Parts.
struct PartByPartIdComparator {
  bool operator()(not_null<Part*> left, not_null<Part*> right) const;
  bool operator()(not_null<Part const*> left,
                  not_null<Part const*> right) const;
};

// Comparator by GUID.  Useful for ensuring a consistent ordering in sets of
// pointers to Vessels.
struct VesselByGUIDComparator {
  bool operator()(not_null<Vessel*> left, not_null<Vessel*> right) const;
  bool operator()(not_null<Vessel const*> left,
                  not_null<Vessel const*> right) const;
};

export template<typename T>
using PartTo = std::map<not_null<Part*>, T, PartByPartIdComparator>;
export using VesselSet = std::set<not_null<Vessel*>, VesselByGUIDComparator>;
export using VesselConstSet =
    std::set<not_null<Vessel const*>, VesselByGUIDComparator>;

}  // namespace ksp_plugin
}  // namespace principia
