// This module contains the other non-SI units listed in the BIPM’s
// SI brochure 8, section 4.1, table 8,
// http://www.bipm.org/en/si/si_brochure/chapter4/table8.html.
export module principia.quantities.bipm;

import principia.quantities.elementary_functions;
import principia.quantities.names;
import principia.quantities.si;

using namespace principia::quantities::elementary_functions;
using namespace principia::quantities::names;
using namespace principia::quantities::si;

export namespace principia::quantities::bipm {

constexpr Pressure Bar                 = 1e5 * Pascal;
constexpr Pressure MillimetreOfMercury = 133.322 * Pascal;
constexpr Length   Ångström            = 1e-10 * Metre;
constexpr Length   NauticalMile        = 1852 * Metre;
constexpr Speed    Knot                = 1 * NauticalMile / Hour;
constexpr Area     Barn                = 1e-28 * Pow<2>(Metre);

}  // namespace principia::quantities::bipm