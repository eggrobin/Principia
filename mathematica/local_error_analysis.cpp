
#include "mathematica/local_error_analysis.hpp"

#include <random>
#include <vector>

#include "astronomy/stabilize_ksp.hpp"
#include "base/bundle.hpp"
#include "base/file.hpp"
#include "mathematica/mathematica.hpp"

namespace principia {
namespace mathematica {

using base::Bundle;
using base::OFStream;
using base::make_not_null_unique;
using base::Status;
using geometry::Position;
using geometry::Vector;
using physics::DegreesOfFreedom;
using physics::MassiveBody;
using quantities::Angle;
using quantities::Cos;
using quantities::Pow;
using quantities::Sin;
using quantities::Sqrt;
using quantities::si::Day;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;

namespace {

constexpr Length fitting_tolerance = 1 * Milli(Metre);

template<typename BitGenerator>
Vector<double, ICRFJ2000Equator> RandomUnitVector(BitGenerator& generator) {
  std::uniform_real_distribution<> longitude_distribution(-π, π);
  std::uniform_real_distribution<> z_distribution(-1, 1);
  double const z = z_distribution(generator);
  Angle const longitude = longitude_distribution(generator) * Radian;
  return Vector<double, ICRFJ2000Equator>({Cos(longitude) * Sqrt(1 - Pow<2>(z)),
                                           Sin(longitude) * Sqrt(1 - Pow<2>(z)),
                                           z});
}

}  // namespace

LocalErrorAnalyser::LocalErrorAnalyser(
    not_null<std::unique_ptr<SolarSystem<ICRFJ2000Equator>>> solar_system,
    FixedStepSizeIntegrator<
        Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const& integrator,
    Time const& step)
    : solar_system_(std::move(solar_system)),
      integrator_(integrator),
      step_(step) {
  // The system might not be defined from Keplerian elements, so we cannot
  // always turn it into a hierarchical system to take its fingerprint.
  // TODO(eggrobin): arbitrary solar system fingerprinting.
  if (solar_system_->names()[0] == "Bop") {
    LOG(INFO) << "All hail retrobop!";
    astronomy::StabilizeKSP(*solar_system_);
  }
}

void LocalErrorAnalyser::WriteLocalErrors(
    std::experimental::filesystem::path const& path,
    FixedStepSizeIntegrator<
        Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
        fine_integrator,
    Time const& fine_step,
    Time const& granularity,
    Time const& begin,
    Time const& end) const {
  auto const reference_ephemeris = solar_system_->MakeEphemeris(
      fitting_tolerance,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(integrator_, step_));
  reference_ephemeris->Prolong(solar_system_->epoch());
  std::vector<std::vector<Length>> errors;
  for (Instant t0 = solar_system_->epoch(),
               t = t0 + granularity;
       t < solar_system_->epoch() + end;
       t0 = t, t += granularity) {
    std::unique_ptr<Ephemeris<ICRFJ2000Equator>> refined_ephemeris;
    reference_ephemeris->Prolong(t);
    if (t >= solar_system_->epoch() + begin) {
      refined_ephemeris =
          ForkEphemeris(*reference_ephemeris, t0, fine_integrator, fine_step);
      refined_ephemeris->Prolong(t);
    }
    LOG_EVERY_N(INFO, 10) << "Prolonged to "
                          << (t - solar_system_->epoch()) / Day << " days.";

    if (t < solar_system_->epoch() + begin) {
      continue;
    }

    errors.emplace_back();
    for (auto const& body_name : solar_system_->names()) {
      int const body_index = solar_system_->index(body_name);
      errors.back().push_back(
          (reference_ephemeris
               ->trajectory(reference_ephemeris->bodies()[body_index])
               ->EvaluatePosition(t) -
           refined_ephemeris
               ->trajectory(refined_ephemeris->bodies()[body_index])
               ->EvaluatePosition(t)).Norm());
    }
  }
  OFStream file(path);
  file << Assign("bodyNames", solar_system_->names());
  file << Assign("errors", ExpressIn(Metre, errors));
}

void LocalErrorAnalyser::WriteEnsembleDiameters(
      std::experimental::filesystem::path const& path,
      Length const& perturbation_norm,
      int const ensemble_size,
      FixedStepSizeIntegrator<
          Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
          fine_integrator,
      Time const& fine_step,
      Time const& granularity,
      Time const& begin,
      Time const& end) const {
  auto const reference_ephemeris = solar_system_->MakeEphemeris(
      fitting_tolerance,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(integrator_, step_));
  reference_ephemeris->Prolong(solar_system_->epoch());
  std::vector<std::vector<Length>> diameters;
  std::vector<Time> times;
  for (Instant t0 = solar_system_->epoch(),
               t = t0 + granularity;
       t < solar_system_->epoch() + end;
       t0 = t, t += granularity) {
    std::vector<not_null<std::unique_ptr<Ephemeris<ICRFJ2000Equator>>>> ensemble;

    Bundle bundle{static_cast<int>(std::thread::hardware_concurrency() - 1)};
    bundle.Add([&reference_ephemeris = *reference_ephemeris, t ]() {
      reference_ephemeris.Prolong(t);
      return Status::OK;
    });
    if (t >= solar_system_->epoch() + begin) {
      ForkEphemerisEnsemble(*reference_ephemeris,
                            t0,
                            perturbation_norm,
                            fine_integrator,
                            fine_step,
                            ensemble_size);
      for (auto const& ephemeris : ensemble) {
        bundle.Add([& ephemeris = *ephemeris, t ]() {
          ephemeris.Prolong(t);
          return Status::OK;
        });
      }
    }
    bundle.Join();
    LOG_EVERY_N(INFO, 10) << "Prolonged to "
                          << (t - solar_system_->epoch()) / Day << " days.";

    if (t < solar_system_->epoch() + begin) {
      continue;
    }

    diameters.emplace_back();
    for (auto const& body_name : solar_system_->names()) {
      Length diameter = 0 * Metre;
      for (int i = 0; i < ensemble.size(); ++i) {
        for (int j = 0; j < i; ++j) {
          int const body_index = solar_system_->index(body_name);
          auto const& trajectory_i =
              *ensemble[i]->trajectory(ensemble[i]->bodies()[body_index]);
          auto const& trajectory_j =
              *ensemble[j]->trajectory(ensemble[j]->bodies()[body_index]);
          diameter = std::max(diameter,
                              (trajectory_i.EvaluatePosition(t) -
                               trajectory_j.EvaluatePosition(t)).Norm());
        }
      }
      diameters.back().push_back(diameter);
    }
    times.push_back(t - solar_system_->epoch());
  }
  OFStream file(path);
  file << Assign("bodyNames", solar_system_->names());
  file << Assign("times", ExpressIn(Second, times));
  file << Assign("diameters", ExpressIn(Metre, diameters));
}

not_null<std::unique_ptr<Ephemeris<ICRFJ2000Equator>>>
LocalErrorAnalyser::ForkEphemeris(
    Ephemeris<ICRFJ2000Equator> const& original,
    Instant const& t,
    FixedStepSizeIntegrator<
        Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const& integrator,
    Time const& step) const {
  std::vector<DegreesOfFreedom<ICRFJ2000Equator>> degrees_of_freedom;
  for (not_null<MassiveBody const*> const body : original.bodies()) {
    degrees_of_freedom.emplace_back(
        original.trajectory(body)->EvaluateDegreesOfFreedom(t));
  }
  return make_not_null_unique<Ephemeris<ICRFJ2000Equator>>(
      solar_system_->MakeAllMassiveBodies(),
      degrees_of_freedom,
      t,
      fitting_tolerance,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(integrator, step));
}

std::vector<not_null<std::unique_ptr<Ephemeris<ICRFJ2000Equator>>>>
LocalErrorAnalyser::ForkEphemerisEnsemble(
    Ephemeris<ICRFJ2000Equator> const& original,
    Instant const& t,
    Length const& perturbation_norm,
    FixedStepSizeIntegrator<
        Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const& integrator,
    Time const& step,
    int const size) const {
  std::mt19937_64 generator;
  std::vector<not_null<std::unique_ptr<Ephemeris<ICRFJ2000Equator>>>> ensemble;
  ensemble.reserve(size);
  for (int i = 0; i < size; ++i) {
    std::vector<DegreesOfFreedom<ICRFJ2000Equator>> degrees_of_freedom;
    degrees_of_freedom.reserve(original.bodies().size());
    for (not_null<MassiveBody const*> const body : original.bodies()) {
      auto const original_degrees_of_freedom =
          original.trajectory(body)->EvaluateDegreesOfFreedom(t);
      auto const perturbed_position =
          original_degrees_of_freedom.position() +
          RandomUnitVector(generator) * perturbation_norm;
      degrees_of_freedom.emplace_back(perturbed_position,
                                      original_degrees_of_freedom.velocity());
    }
    ensemble.emplace_back(make_not_null_unique<Ephemeris<ICRFJ2000Equator>>(
        solar_system_->MakeAllMassiveBodies(),
        degrees_of_freedom,
        t,
        fitting_tolerance,
        Ephemeris<ICRFJ2000Equator>::FixedStepParameters(integrator, step)));
  }
  return ensemble;
}

}  // namespace mathematica
}  // namespace principia
