
#pragma once

#include <vector>

#include "integrators/integrators.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"

namespace principia {
namespace mathematica {

using astronomy::ICRFJ2000Equator;
using base::not_null;
using geometry::Instant;
using integrators::FixedStepSizeIntegrator;
using physics::Ephemeris;
using physics::SolarSystem;
using quantities::Length;
using quantities::Time;

// A utility to compute the local errors in the numerical integration of a
// |solar_system| with a given |integrator| and |step|.
class LocalErrorAnalyser {
 public:
  LocalErrorAnalyser(
      not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>> solar_system,
      FixedStepSizeIntegrator<
          Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
          integrator,
      Time const& step);

  // Computes the error over |granularity| between the main integration and a
  // fine integration forked off the main one, for |duration| from the solar
  // system epoch.  Writes the errors to a file with the given |path|.
  void WriteLocalErrors(
      std::experimental::filesystem::path const& path,
      FixedStepSizeIntegrator<
          Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
          fine_integrator,
      Time const& fine_step,
      Time const& granularity,
      Time const& duration) const;

  void WriteEnsembleDiameters(
      std::experimental::filesystem::path const& path,
      Length const& perturbation_norm,
      int ensemble_size,
      FixedStepSizeIntegrator<
          Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
          fine_integrator,
      Time const& fine_step,
      Time const& granularity,
      Time const& duration) const;

 private:
  not_null<std::unique_ptr<Ephemeris<ICRFJ2000Equator>>> ForkEphemeris(
      Ephemeris<ICRFJ2000Equator> const& original,
      Instant const& t,
      FixedStepSizeIntegrator<
          Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
          integrator,
      Time const& step) const;

  std::vector<not_null<std::unique_ptr<Ephemeris<ICRFJ2000Equator>>>>
  ForkEphemerisEnsemble(
      Ephemeris<ICRFJ2000Equator> const& original,
      Instant const& t,
      Length const& perturbation_norm,
      FixedStepSizeIntegrator<
          Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
          integrator,
      Time const& step,
      int size) const;

  not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>> const solar_system_;
  FixedStepSizeIntegrator<
      Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const& integrator_;
  Time const step_;
};

}  // namespace mathematica
}  // namespace principia
