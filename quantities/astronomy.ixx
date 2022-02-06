// This module contains units commonly used in astronomy.
export module principia.quantities.astronomy;

import "quantities/numbers.hpp";

import principia.quantities.constants;
import principia.quantities.elementary_functions;
import principia.quantities.names;
import principia.quantities.si;

using namespace principia::quantities::constants;
using namespace principia::quantities::elementary_functions;
using namespace principia::quantities::names;
using namespace principia::quantities::si;

export namespace principia::quantities::astronomy {

// Résolution B2 "Re-définition de l’unité astronomique de longueur" adopted
// at the XXVIIIth General Assembly of the IAU in 2012.
constexpr Length AstronomicalUnit = 149597870700 * Metre;

// See note 4 of Résolution B2 "Sur la recommandation du "point zéro" des
// échelles de magnitude bolométrique absolue et apparente" adopted at the
// XXIXth General Assembly of the IAU in 2015.
constexpr Length Parsec = 648000 / π * AstronomicalUnit;

// System of nominal solar and planetary conversion constants, Résolution B3
// "Sur les valeurs  recommandées de constantes de conversion pour une sélection
// de propriétés solaires et planétaires" adopted at the XXIXth General Assembly
// of the IAU in 2015.

// Solar conversion constants.
constexpr Length      SolarRadius               = 6.957e8 * Metre;
constexpr Irradiance  TotalSolarIrradiance      = 1361 * (Watt /
                                                          Pow<2>(Metre));
constexpr Power       SolarLuminosity           = 3.828e26 * Watt;
constexpr Temperature SolarEffectiveTemperature = 5772 * Kelvin;
constexpr GravitationalParameter SolarGravitationalParameter =
    1.327'124'4e20 * (Pow<3>(Metre) / Pow<2>(Second));

// Planetary conversion constants.
// “If equatorial vs. polar radius is not explicitly specified, it should be
// understood that nominal terrestrial [or jovian] radius refers specifically to
// [the nominal equatorial radius], following common usage.”
constexpr Length TerrestrialEquatorialRadius = 6.3781e6 * Metre;
constexpr Length TerrestrialPolarRadius      = 6.3568e6 * Metre;
constexpr Length JovianEquatorialRadius      = 7.1492e7 * Metre;
constexpr Length JovianPolarRadius           = 6.6854e7 * Metre;
constexpr GravitationalParameter TerrestrialGravitationalParameter =
    3.986'004e14 * (Pow<3>(Metre) / Pow<2>(Second));
constexpr GravitationalParameter JovianGravitationalParameter      =
    1.266'865'3e17 * (Pow<3>(Metre) / Pow<2>(Second));

constexpr Time   JulianYear = 365.25 * Day;
constexpr Length LightYear  = SpeedOfLight * JulianYear;

}  // namespace principia::quantities::astronomy
