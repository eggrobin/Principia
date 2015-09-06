﻿#pragma once

#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/manœuvre.hpp"
#include "ksp_plugin/plugin.hpp"

namespace principia {
namespace ksp_plugin {

// TODO(egg): that constructor is going to get annoying, set everything with
// single-use mutators.
extern "C" DLLEXPORT
Manœuvre* principia__NewManœuvreIspByWeight(
    double thrust,
    double initial_mass,
    double specific_impulse_by_weight,
    double right_ascension,
    double declination);

extern "C" DLLEXPORT
double principia__thrust(Manœuvre const* manœuvre);
extern "C" DLLEXPORT
double principia__initial_mass(Manœuvre const* manœuvre);
extern "C" DLLEXPORT
double principia__specific_impulse_by_weight(
    Manœuvre const* manœuvre);

extern "C" DLLEXPORT
XYZ principia__Δv(Manœuvre const* manœuvre);

extern "C" DLLEXPORT
void principia__set_duration(Manœuvre* manœuvre, double duration);
extern "C" DLLEXPORT
void principia__set_Δv(Manœuvre const* manœuvre, double Δv);

extern "C" DLLEXPORT
double principia__initial_time(Manœuvre const* manœuvre);
extern "C" DLLEXPORT
void principia__set_initial_time(Manœuvre* manœuvre, double initial_time);
extern "C" DLLEXPORT
double principia__time_of_half_Δv(Manœuvre const* manœuvre);
extern "C" DLLEXPORT
void principia__set_time_of_half_Δv(Manœuvre* manœuvre, double time_of_half_Δv);

}  // namespace ksp_plugin
}  // namespace principia
