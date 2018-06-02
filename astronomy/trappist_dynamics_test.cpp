
#include <random>

#include "astronomy/frames.hpp"
#include "base/bundle.hpp"
#include "base/file.hpp"
#include "base/not_null.hpp"
#include "base/status.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/sign.hpp"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/root_finders.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::Bundle;
using base::not_null;
using base::OFStream;
using base::Status;
using geometry::Instant;
using geometry::Position;
using geometry::Sign;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using numerics::Bisect;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::RelativeDegreesOfFreedom;
using physics::SolarSystem;
using quantities::Angle;
using quantities::Derivative;
using quantities::Difference;
using quantities::Exponentiation;
using quantities::Pow;
using quantities::Square;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;

namespace astronomy {

struct PlanetParameters {
  constexpr static int count = 4;
  double eccentricity;
  Time period;
  Angle argument_of_periapsis;
  Angle mean_anomaly;
};

using SystemParameters = std::array<PlanetParameters, 7>;

using Population = std::vector<SystemParameters>;

using Calculator = std::function<double(SystemParameters const&)>;

PlanetParameters operator+(PlanetParameters const& left,
                           PlanetParameters const& right) {
  PlanetParameters result = left;
  result.eccentricity += right.eccentricity;
  result.period += right.period;
  result.argument_of_periapsis += right.argument_of_periapsis;
  result.mean_anomaly += right.mean_anomaly;
  return result;
}

PlanetParameters operator-(PlanetParameters const& left,
                           PlanetParameters const& right) {
  PlanetParameters result = left;
  result.eccentricity -= right.eccentricity;
  result.period -= right.period;
  result.argument_of_periapsis -= right.argument_of_periapsis;
  result.mean_anomaly -= right.mean_anomaly;
  return result;
}

PlanetParameters operator*(double const left, PlanetParameters const& right) {
  PlanetParameters result = right;
  result.eccentricity *= left;
  result.period *= left;
  result.argument_of_periapsis *= left;
  result.mean_anomaly *= left;
  return result;
}

SystemParameters operator+(SystemParameters const& left,
                           SystemParameters const& right) {
  SystemParameters result = left;
  for (int i = 0; i < result.size(); ++i) {
    result[i] = result[i] + right[i];
  }
  return result;
}

SystemParameters operator-(SystemParameters const& left,
                           SystemParameters const& right) {
  SystemParameters result = left;
  for (int i = 0; i < result.size(); ++i) {
    result[i] = result[i] - right[i];
  }
  return result;
}

SystemParameters operator*(double const left, SystemParameters const& right) {
  SystemParameters result = right;
  for (int i = 0; i < result.size(); ++i) {
    result[i] = left * result[i];
  }
  return result;
}

KeplerianElements<Trappist> MakeKeplerianElements(
    KeplerianElements<Trappist> const& blueprint,
    PlanetParameters const& parameters) {
  KeplerianElements<Trappist> elements = blueprint;
  elements.asymptotic_true_anomaly = std::nullopt;
  elements.turning_angle = std::nullopt;
  elements.semimajor_axis = std::nullopt;
  elements.specific_energy = std::nullopt;
  elements.characteristic_energy = std::nullopt;
  elements.mean_motion = std::nullopt;
  elements.hyperbolic_mean_motion = std::nullopt;
  elements.hyperbolic_excess_velocity = std::nullopt;
  elements.semiminor_axis = std::nullopt;
  elements.impact_parameter = std::nullopt;
  elements.semilatus_rectum = std::nullopt;
  elements.specific_angular_momentum = std::nullopt;
  elements.periapsis_distance = std::nullopt;
  elements.apoapsis_distance = std::nullopt;
  elements.longitude_of_periapsis = std::nullopt;
  elements.true_anomaly = std::nullopt;
  elements.hyperbolic_mean_anomaly = std::nullopt;
  *elements.argument_of_periapsis = parameters.argument_of_periapsis;
  *elements.mean_anomaly = parameters.mean_anomaly;
  *elements.period = parameters.period;
  *elements.eccentricity = parameters.eccentricity;
  return elements;
}

PlanetParameters MakePlanetParameters(
    KeplerianElements<Trappist> const& elements) {
  PlanetParameters result;
  result.argument_of_periapsis = *elements.argument_of_periapsis;
  result.eccentricity = *elements.eccentricity;
  result.mean_anomaly = *elements.mean_anomaly;
  result.period = *elements.period;
  return result;
}

std::vector<double> EvaluatePopulation(
    Population const& population,
    Calculator const& calculate_log_pdf) {
  std::vector<double> log_pdf(population.size());
  Bundle bundle(4);
  for (int i = 0; i < population.size(); ++i) {
    auto const& parameters = population[i];
    bundle.Add([&calculate_log_pdf, i, &log_pdf, &parameters]() {
      log_pdf[i] = calculate_log_pdf(parameters);
      return Status::OK;
    });
  }
  bundle.Join();
  return log_pdf;
}

Population GenerateTrialStatesDEMCMC(Population const& population,
                                     double const γ,
                                     double const ε,
                                     std::mt19937_64& engine) {
  Population trial(population.size());
  std::uniform_int_distribution<> j_distribution(0, population.size() - 2);
  std::uniform_int_distribution<> k_distribution(0, population.size() - 3);
  std::normal_distribution<> perturbation_distribution;
  for (int i = 0; i < population.size(); ++i) {
    // Choose head (k) and tail (j) for perturbation vector.
    int j = j_distribution(engine);
    if (j >= i) {
      ++j;
    }
    int k = k_distribution(engine);
    if (k >= i) {
      ++k;
    }
    if (k >= j) {
      ++k;
    }

    // Choose scale factor.
    double const scale = (1.0 + ε * perturbation_distribution(engine)) * γ;

    trial[i] = population[i] + scale * (population[k] - population[j]);
  }
  return trial;
}

void RunDEMCMC(Population& population,
               int const number_of_generations,
               int const number_of_generations_between_kicks,
               int const number_of_burn_in_generations,
               double const ε,
               Calculator const& calculate_log_pdf) {
  CHECK_LE(1, number_of_generations);
  CHECK_LT(std::tuple_size<SystemParameters>::value * PlanetParameters::count,
           population.size());
  std::mt19937_64 engine;
  std::uniform_real_distribution<> distribution(0.0, 1.0);

  static double best_log_pdf = std::numeric_limits<double>::lowest();
  std::vector<Population> ancestry;
  std::vector<std::vector<double>> ancestry_log_pdf;
  std::vector<int> accepts_chain(population.size(), 0);
  std::vector<int> rejects_chain(population.size(), 0);
  std::vector<int> accepts_generation(number_of_generations, 0);
  std::vector<int> rejects_generation(number_of_generations, 0);

  auto log_pdf = EvaluatePopulation(population, calculate_log_pdf);

  // Loop over generations.
  for (int generation = 0; generation < number_of_generations; ++generation) {
    LOG_IF(ERROR, generation % 50 == 0) << "Generation: " << generation;

    // Every 10th generation try full-size steps.
    double const
        γ = generation < number_of_burn_in_generations
                ? 0.01
                : generation % number_of_generations_between_kicks == 0
                      ? 1.0
                      : 2.38 /
                            Sqrt(2 * std::tuple_size<SystemParameters>::value *
                                 PlanetParameters::count);

    // Evaluate model for each set of trial parameters.
    auto const trial = GenerateTrialStatesDEMCMC(population, γ, ε, engine);
    auto const log_pdf_trial = EvaluatePopulation(trial, calculate_log_pdf);

    // For each member of population.
    for (int i = 0; i < population.size(); ++i) {
      double const log_pdf_ratio = log_pdf_trial[i] - log_pdf[i];
      if (log_pdf_ratio > 0.0 ||
          log_pdf_ratio > std::log(distribution(engine))) {
        population[i] = trial[i];
        log_pdf[i] = log_pdf_trial[i];
        ++accepts_chain[i];
        ++accepts_generation[generation];
      } else {
        ++rejects_chain[i];
        ++rejects_generation[generation];
      }
    }

    // Traces.
    best_log_pdf = std::max(best_log_pdf,
                            *std::max_element(log_pdf.begin(), log_pdf.end()));
    LOG(ERROR) << "Min: " << *std::min_element(log_pdf.begin(), log_pdf.end())
               << " Max: " << *std::max_element(log_pdf.begin(), log_pdf.end())
               << " Best: " << best_log_pdf;

    // Record results.
    ancestry.push_back(population);
    ancestry_log_pdf.push_back(log_pdf);
  }

  for (int generation = 0; generation < number_of_generations; ++generation) {
    LOG(ERROR) << generation << " Accepts: " << accepts_generation[generation]
               << " Rejects: " << rejects_generation[generation];
  }
}

#if 0
void Genome::SubStepOptimization(
    std::function<double(Genome const&)> const& χ²) {
  Bundle bundle(4);
  double χ²₁;
  bundle.Add([this, &χ², &χ²₁]() {
    χ²₁ = χ²(*this);
    return Status::OK;
  });
  std::vector<Angle> M₂s;
  M₂s.resize(elements_.size());
  std::vector<Angle> M₃s;
  M₃s.resize(elements_.size());
  std::vector<double> χ²₂s;
  χ²₂s.resize(elements_.size());
  std::vector<double> χ²₃s;
  χ²₃s.resize(elements_.size());
  for (int i = 0; i < elements_.size(); ++i) {
    Genome perturbed_genome = *this;
    auto& perturbed_element = perturbed_genome.elements_[i];
    *perturbed_element.mean_anomaly += 0.1 * Degree;
    M₂s[i] = *perturbed_element.mean_anomaly;
    bundle.Add([perturbed_genome, i, &χ², &χ²₂s]() {
      χ²₂s [i] = χ²(perturbed_genome);
      return Status::OK;
    });
    *perturbed_element.mean_anomaly -= 0.2 * Degree;
    M₃s[i] = *perturbed_element.mean_anomaly;
      bundle.Add([perturbed_genome, i, &χ², &χ²₃s]() {
      χ²₃s [i] = χ²(perturbed_genome);
      return Status::OK;
    });
  }
  bundle.Join();

  for (int i = 0; i < elements_.size(); ++i) {
    double const χ²₂ = χ²₂s [i];
    double const χ²₃ = χ²₃s [i];
    Angle const M₁ = *elements_[i].mean_anomaly;
    Angle const M₂ = M₂s[i];
    Angle const M₃ = M₃s[i];
    Derivative<double, Angle> b₁ = (χ²₂ - χ²₁) / (M₂ - M₁);
    Derivative<double, Angle, 2> b₂ =
        ((χ²₃ - χ²₁) / (M₃ - M₁) - b₁) / (M₃ - M₂);
    elements_[i].mean_anomaly = (M₁ + M₂) / 2 - b₁ / (2 * b₂);
  }
}
#endif

// TODO(phl): Literals are broken in 15.8.0 Preview 1.0 and are off by an
// integral number of days.  Use this function as a stopgap measure and switch
// to literals once MSFT have fixed their bugs.
constexpr Instant JD(double const jd) {
  return Instant{} + (jd - 2451545.0) * Day;
}

double ShortDays(Instant const& time) {
    return (time - JD(2450000.0)) / Day;
}

using Transits = std::vector<Instant>;
using TransitsByPlanet = std::map<std::string, Transits>;

template<typename Value>
struct Measured {
  Value estimated_value;
  Difference<Value> standard_uncertainty;
};

using MeasuredTransits = std::vector<Measured<Instant>>;
using MeasuredTransitsByPlanet = std::map<std::string, MeasuredTransits>;

MeasuredTransitsByPlanet const observations = {
    {"Trappist-1b",
     {{JD(2457322.51531), 0.00071 * Day}, {JD(2457325.53910), 0.00100 * Day},
      {JD(2457328.55860), 0.00130 * Day}, {JD(2457331.58160), 0.00100 * Day},
      {JD(2457334.60480), 0.00017 * Day}, {JD(2457337.62644), 0.00092 * Day},
      {JD(2457340.64820), 0.00140 * Day}, {JD(2457345.18028), 0.00080 * Day},
      {JD(2457361.79945), 0.00028 * Day}, {JD(2457364.82173), 0.00077 * Day},
      {JD(2457440.36492), 0.00020 * Day}, {JD(2457452.45228), 0.00014 * Day},
      {JD(2457463.02847), 0.00019 * Day}, {JD(2457509.86460), 0.00210 * Day},
      {JD(2457512.88731), 0.00029 * Day}, {JD(2457568.78880), 0.00100 * Day},
      {JD(2457586.91824), 0.00064 * Day}, {JD(2457589.93922), 0.00092 * Day},
      {JD(2457599.00640), 0.00021 * Day}, {JD(2457602.02805), 0.00071 * Day},
      {JD(2457612.60595), 0.00085 * Day}, {JD(2457615.62710), 0.00160 * Day},
      {JD(2457624.69094), 0.00066 * Day}, {JD(2457645.84400), 0.00110 * Day},
      {JD(2457651.88743), 0.00022 * Day}, {JD(2457653.39809), 0.00026 * Day},
      {JD(2457654.90908), 0.00084 * Day}, {JD(2457656.41900), 0.00029 * Day},
      {JD(2457657.93129), 0.00020 * Day}, {JD(2457659.44144), 0.00017 * Day},
      {JD(2457660.95205), 0.00035 * Day}, {JD(2457662.46358), 0.00020 * Day},
      {JD(2457663.97492), 0.00070 * Day}, {JD(2457665.48509), 0.00017 * Day},
      {JD(2457666.99567), 0.00025 * Day}, {JD(2457668.50668), 0.00030 * Day},
      {JD(2457670.01766), 0.00034 * Day}, {JD(2457671.52876), 0.00033 * Day},
      {JD(2457721.38747), 0.00035 * Day}, {JD(2457739.51770), 0.00059 * Day},
      {JD(2457741.02787), 0.00055 * Day}, {JD(2457742.53918), 0.00058 * Day},
      {JD(2457744.05089), 0.00061 * Day}, {JD(2457745.56164), 0.00072 * Day},
      {JD(2457747.07208), 0.00085 * Day}, {JD(2457748.58446), 0.00087 * Day},
      {JD(2457750.09387), 0.00089 * Day}, {JD(2457751.60535), 0.00082 * Day},
      {JD(2457753.11623), 0.00075 * Day}, {JD(2457754.62804), 0.00077 * Day},
      {JD(2457756.13856), 0.00060 * Day}, {JD(2457757.64840), 0.00089 * Day},
      {JD(2457759.15953), 0.00073 * Day}, {JD(2457760.67112), 0.00082 * Day},
      {JD(2457762.18120), 0.00073 * Day}, {JD(2457763.69221), 0.00071 * Day},
      {JD(2457765.20298), 0.00077 * Day}, {JD(2457766.71479), 0.00055 * Day},
      {JD(2457768.22514), 0.00103 * Day}, {JD(2457769.73704), 0.00064 * Day},
      {JD(2457771.24778), 0.00091 * Day}, {JD(2457772.75738), 0.00075 * Day},
      {JD(2457774.26841), 0.00080 * Day}, {JD(2457775.77995), 0.00058 * Day},
      {JD(2457777.28899), 0.00099 * Day}, {JD(2457778.80118), 0.00062 * Day},
      {JD(2457780.31297), 0.00068 * Day}, {JD(2457781.82231), 0.00145 * Day},
      {JD(2457783.33410), 0.00071 * Day}, {JD(2457784.84372), 0.00068 * Day},
      {JD(2457792.39979), 0.00110 * Day}, {JD(2457793.90955), 0.00064 * Day},
      {JD(2457795.41987), 0.00058 * Day}, {JD(2457796.93134), 0.00065 * Day},
      {JD(2457798.44211), 0.00061 * Day}, {JD(2457799.95320), 0.00083 * Day},
      {JD(2457801.46314), 0.00127 * Day}, {JD(2457802.97557), 0.00016 * Day},
      {JD(2457804.48638), 0.00053 * Day}, {JD(2457805.99697), 0.00016 * Day},
      {JD(2457807.50731), 0.00017 * Day}, {JD(2457809.01822), 0.00017 * Day},
      {JD(2457810.52781), 0.00110 * Day}, {JD(2457812.04038), 0.00020 * Day},
      {JD(2457813.55121), 0.00014 * Day}, {JD(2457815.06275), 0.00017 * Day},
      {JD(2457816.57335), 0.00011 * Day}, {JD(2457818.08382), 0.00015 * Day},
      {JD(2457819.59478), 0.00017 * Day}, {JD(2457821.10550), 0.00020 * Day},
      {JD(2457824.12730), 0.00018 * Day}, {JD(2457825.63813), 0.00018 * Day},
      {JD(2457827.14995), 0.00012 * Day}, {JD(2457828.66042), 0.00024 * Day},
      {JD(2457830.17087), 0.00021 * Day}, {JD(2457833.19257), 0.00018 * Day},
      {JD(2457834.70398), 0.00016 * Day}, {JD(2457836.21440), 0.00017 * Day},
      {JD(2457837.72526), 0.00014 * Day}, {JD(2457839.23669), 0.00017 * Day},
      {JD(2457917.80060), 0.00110 * Day}, {JD(2457923.84629), 0.00045 * Day},
      {JD(2457935.93288), 0.00023 * Day}, {JD(2457952.55450), 0.00110 * Day},
      {JD(2457955.57554), 0.00069 * Day}, {JD(2457967.66254), 0.00050 * Day},
      {JD(2457973.70596), 0.00040 * Day}}},
    {"Trappist-1c",
     {{JD(2457282.80570), 0.00140 * Day}, {JD(2457333.66400), 0.00090 * Day},
      {JD(2457362.72605), 0.00038 * Day}, {JD(2457367.57051), 0.00033 * Day},
      {JD(2457384.52320), 0.00130 * Day}, {JD(2457452.33470), 0.00015 * Day},
      {JD(2457454.75672), 0.00066 * Day}, {JD(2457512.88094), 0.00009 * Day},
      {JD(2457546.78587), 0.00075 * Day}, {JD(2457551.62888), 0.00066 * Day},
      {JD(2457580.69137), 0.00031 * Day}, {JD(2457585.53577), 0.00250 * Day},
      {JD(2457587.95622), 0.00054 * Day}, {JD(2457600.06684), 0.00036 * Day},
      {JD(2457604.90975), 0.00063 * Day}, {JD(2457609.75461), 0.00072 * Day},
      {JD(2457614.59710), 0.00130 * Day}, {JD(2457626.70610), 0.00110 * Day},
      {JD(2457631.55024), 0.00056 * Day}, {JD(2457638.81518), 0.00048 * Day},
      {JD(2457650.92395), 0.00023 * Day}, {JD(2457653.34553), 0.00024 * Day},
      {JD(2457655.76785), 0.00043 * Day}, {JD(2457658.18963), 0.00024 * Day},
      {JD(2457660.61168), 0.00051 * Day}, {JD(2457663.03292), 0.00028 * Day},
      {JD(2457665.45519), 0.00025 * Day}, {JD(2457667.87729), 0.00031 * Day},
      {JD(2457670.29869), 0.00035 * Day}, {JD(2457672.71944), 0.00081 * Day},
      {JD(2457711.46778), 0.00064 * Day}, {JD(2457723.57663), 0.00050 * Day},
      {JD(2457740.53361), 0.00088 * Day}, {JD(2457742.95276), 0.00115 * Day},
      {JD(2457745.37429), 0.00063 * Day}, {JD(2457747.79699), 0.00056 * Day},
      {JD(2457750.21773), 0.00096 * Day}, {JD(2457752.64166), 0.00093 * Day},
      {JD(2457755.05877), 0.00165 * Day}, {JD(2457757.48313), 0.00066 * Day},
      {JD(2457759.90281), 0.00058 * Day}, {JD(2457762.32806), 0.00081 * Day},
      {JD(2457764.74831), 0.00072 * Day}, {JD(2457767.16994), 0.00125 * Day},
      {JD(2457769.59209), 0.00081 * Day}, {JD(2457772.01483), 0.00100 * Day},
      {JD(2457774.43458), 0.00081 * Day}, {JD(2457776.85815), 0.00102 * Day},
      {JD(2457779.27911), 0.00089 * Day}, {JD(2457781.70095), 0.00072 * Day},
      {JD(2457784.12338), 0.00054 * Day}, {JD(2457791.38801), 0.00064 * Day},
      {JD(2457793.81141), 0.00079 * Day}, {JD(2457796.23153), 0.00052 * Day},
      {JD(2457798.65366), 0.00082 * Day}, {JD(2457801.07631), 0.00084 * Day},
      {JD(2457803.49747), 0.00020 * Day}, {JD(2457805.91882), 0.00017 * Day},
      {JD(2457808.34123), 0.00023 * Day}, {JD(2457810.76273), 0.00019 * Day},
      {JD(2457813.18456), 0.00024 * Day}, {JD(2457815.60583), 0.00017 * Day},
      {JD(2457818.02821), 0.00020 * Day}, {JD(2457820.45019), 0.00022 * Day},
      {JD(2457822.87188), 0.00021 * Day}, {JD(2457825.29388), 0.00022 * Day},
      {JD(2457827.71513), 0.00022 * Day}, {JD(2457830.13713), 0.00026 * Day},
      {JD(2457832.55888), 0.00015 * Day}, {JD(2457834.98120), 0.00025 * Day},
      {JD(2457837.40280), 0.00017 * Day}, {JD(2457839.82415), 0.00031 * Day}}},
    {"Trappist-1d",
     {{JD(2457560.79730), 0.00230 * Day}, {JD(2457625.59779), 0.00078 * Day},
      {JD(2457641.79360), 0.00290 * Day}, {JD(2457645.84360), 0.00210 * Day},
      {JD(2457653.94261), 0.00051 * Day}, {JD(2457657.99220), 0.00063 * Day},
      {JD(2457662.04284), 0.00051 * Day}, {JD(2457666.09140), 0.00160 * Day},
      {JD(2457670.14198), 0.00066 * Day}, {JD(2457726.83975), 0.00029 * Day},
      {JD(2457738.99169), 0.00160 * Day}, {JD(2457743.03953), 0.00180 * Day},
      {JD(2457747.08985), 0.00145 * Day}, {JD(2457751.14022), 0.00195 * Day},
      {JD(2457755.18894), 0.00155 * Day}, {JD(2457759.24638), 0.00225 * Day},
      {JD(2457763.28895), 0.00150 * Day}, {JD(2457767.33866), 0.00190 * Day},
      {JD(2457771.39077), 0.00260 * Day}, {JD(2457775.44026), 0.00125 * Day},
      {JD(2457779.48843), 0.00190 * Day}, {JD(2457783.54023), 0.00240 * Day},
      {JD(2457791.64083), 0.00135 * Day}, {JD(2457803.79083), 0.00049 * Day},
      {JD(2457807.84032), 0.00030 * Day}, {JD(2457811.89116), 0.00050 * Day},
      {JD(2457815.94064), 0.00030 * Day}, {JD(2457819.99050), 0.00050 * Day},
      {JD(2457824.04185), 0.00067 * Day}, {JD(2457828.09082), 0.00043 * Day},
      {JD(2457832.14036), 0.00037 * Day}, {JD(2457836.19171), 0.00042 * Day},
      {JD(2457961.73760), 0.00130 * Day}, {JD(2457969.83708), 0.00068 * Day},
      {JD(2457973.88590), 0.00066 * Day}}},
    {"Trappist-1e",
     {{JD(2457312.71300), 0.00270 * Day}, {JD(2457367.59683), 0.00037 * Day},
      {JD(2457611.57620), 0.00310 * Day}, {JD(2457623.77950), 0.00100 * Day},
      {JD(2457654.27862), 0.00049 * Day}, {JD(2457660.38016), 0.00078 * Day},
      {JD(2457666.48030), 0.00180 * Day}, {JD(2457672.57930), 0.00260 * Day},
      {JD(2457721.37514), 0.00099 * Day}, {JD(2457733.57300), 0.00140 * Day},
      {JD(2457739.67085), 0.00135 * Day}, {JD(2457745.77160), 0.00120 * Day},
      {JD(2457751.87007), 0.00034 * Day}, {JD(2457757.96712), 0.00160 * Day},
      {JD(2457764.06700), 0.00240 * Day}, {JD(2457770.17109), 0.00215 * Day},
      {JD(2457776.26378), 0.00160 * Day}, {JD(2457782.36226), 0.00175 * Day},
      {JD(2457794.56159), 0.00160 * Day}, {JD(2457800.66354), 0.00170 * Day},
      {JD(2457806.75758), 0.00041 * Day}, {JD(2457812.85701), 0.00034 * Day},
      {JD(2457818.95510), 0.00030 * Day}, {JD(2457825.05308), 0.00035 * Day},
      {JD(2457831.15206), 0.00027 * Day}, {JD(2457837.24980), 0.00025 * Day},
      {JD(2457934.83095), 0.00050 * Day}, {JD(2457940.92995), 0.00086 * Day}}},
    {"Trappist-1f",
     {{JD(2457321.52520), 0.00200 * Day}, {JD(2457367.57629), 0.00044 * Day},
      {JD(2457634.57809), 0.00061 * Day}, {JD(2457652.98579), 0.00032 * Day},
      {JD(2457662.18747), 0.00040 * Day}, {JD(2457671.39279), 0.00072 * Day},
      {JD(2457717.41541), 0.00091 * Day}, {JD(2457726.61960), 0.00026 * Day},
      {JD(2457745.03116), 0.00135 * Day}, {JD(2457754.23380), 0.00155 * Day},
      {JD(2457763.44338), 0.00024 * Day}, {JD(2457772.64752), 0.00160 * Day},
      {JD(2457781.85142), 0.00180 * Day}, {JD(2457800.27307), 0.00140 * Day},
      {JD(2457809.47554), 0.00027 * Day}, {JD(2457818.68271), 0.00032 * Day},
      {JD(2457827.88669), 0.00030 * Day}, {JD(2457837.10322), 0.00032 * Day},
      {JD(2457956.80549), 0.00054 * Day}}},
    {"Trappist-1g",
     {{JD(2457294.78600), 0.00390 * Day}, {JD(2457356.53410), 0.00200 * Day},
      {JD(2457615.92400), 0.00170 * Day}, {JD(2457640.63730), 0.00100 * Day},
      {JD(2457652.99481), 0.00030 * Day}, {JD(2457665.35151), 0.00028 * Day},
      {JD(2457739.48441), 0.00115 * Day}, {JD(2457751.83993), 0.00017 * Day},
      {JD(2457764.19098), 0.00155 * Day}, {JD(2457776.54900), 0.00110 * Day},
      {JD(2457801.25000), 0.00093 * Day}, {JD(2457813.60684), 0.00023 * Day},
      {JD(2457825.96112), 0.00020 * Day}, {JD(2457838.30655), 0.00028 * Day},
      {JD(2457924.77090), 0.00140 * Day}, {JD(2457961.82621), 0.00068 * Day}}},
    {"Trappist-1h",
     {{JD(2457662.55467), 0.00054 * Day}, {JD(2457756.38740), 0.00130 * Day},
      {JD(2457775.15390), 0.00160 * Day}, {JD(2457793.92300), 0.00250 * Day},
      {JD(2457812.69870), 0.00450 * Day}, {JD(2457831.46625), 0.00047 * Day},
      {JD(2457962.86271), 0.00083 * Day}}}};

double Transitsχ²(MeasuredTransitsByPlanet const& observations,
                  TransitsByPlanet const& computations,
                  bool const verbose) {
  double sum_of_squared_errors = 0;
  Time max_error;
  for (auto const& pair : observations) {
    auto const& name = pair.first;
    auto const& observed_transits = pair.second;
    auto const& computed_transits = computations.at(name);
    if (computed_transits.empty()) {
      return std::numeric_limits<double>::infinity();
    }
    for (auto const& observed_transit : observed_transits) {
      auto const next_computed_transit =
          std::lower_bound(computed_transits.begin(),
                           computed_transits.end(),
                           observed_transit.estimated_value);
      Time error;
      if (next_computed_transit == computed_transits.begin()) {
        error = *next_computed_transit - observed_transit.estimated_value;
      } else if (next_computed_transit == computed_transits.end()) {
        error = observed_transit.estimated_value - computed_transits.back();
      } else {
        error =
            std::min(*next_computed_transit - observed_transit.estimated_value,
                     observed_transit.estimated_value -
                         *std::prev(next_computed_transit));
      }
      CHECK_LE(0.0 * Second, error);
      LOG_IF(ERROR, verbose) << name << ": " << error;
      if (error > max_error) {
        max_error = error;
        LOG_IF(ERROR, verbose)
            << name << ": " << ShortDays(*std::prev(next_computed_transit))
            << " " << ShortDays(observed_transit.estimated_value) << " "
            << ShortDays(*next_computed_transit) << " " << error;
      }
      sum_of_squared_errors +=
          Pow<2>(error / observed_transit.standard_uncertainty);
    }
  }
  return sum_of_squared_errors;
}

class TrappistDynamicsTest : public ::testing::Test {
 protected:
  TrappistDynamicsTest()
      : system_(SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
                SOLUTION_DIR / "astronomy" /
                    "trappist_initial_state_jd_2457010_000000000.proto.txt"),
        ephemeris_(system_.MakeEphemeris(
            /*fitting_tolerance=*/5 * Milli(Metre),
            Ephemeris<Trappist>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                                   Position<Trappist>>(),
                /*step=*/0.07 * Day))) {}

  static Transits ComputeTransits(Ephemeris<Trappist> const& ephemeris,
                                  not_null<MassiveBody const*> const star,
                                  not_null<MassiveBody const*> const planet) {
    Transits transits;
    auto const& star_trajectory = ephemeris.trajectory(star);

    std::optional<Instant> last_t;
    std::optional<Sign> last_xy_displacement_derivative_sign;
      auto const& planet_trajectory = ephemeris.trajectory(planet);
    for (Instant t = ephemeris.t_min();
          t < ephemeris.t_max();
          t += 2 * Hour) {
      RelativeDegreesOfFreedom<Trappist> const relative_dof =
          planet_trajectory->EvaluateDegreesOfFreedom(t) -
          star_trajectory->EvaluateDegreesOfFreedom(t);

      auto const xy_displacement_derivative =
          [&planet_trajectory, &star_trajectory](Instant const& t) {
            RelativeDegreesOfFreedom<Trappist> const relative_dof =
                planet_trajectory->EvaluateDegreesOfFreedom(t) -
                star_trajectory->EvaluateDegreesOfFreedom(t);
            // TODO(phl): Why don't we have projections?
            auto xy_displacement =
                relative_dof.displacement().coordinates();
            xy_displacement.z = 0.0 * Metre;
            auto xy_velocity = relative_dof.velocity().coordinates();
            xy_velocity.z = 0.0 * Metre / Second;
            return Dot(xy_displacement, xy_velocity);
          };

      Sign const xy_displacement_derivative_sign(
          xy_displacement_derivative(t));
      if (relative_dof.displacement().coordinates().z > 0.0 * Metre &&
          last_t &&
          xy_displacement_derivative_sign == Sign(1) &&
          last_xy_displacement_derivative_sign == Sign(-1)) {
        Instant const transit =
            Bisect(xy_displacement_derivative, *last_t, t);
        transits.push_back(transit);
      }
      last_t = t;
      last_xy_displacement_derivative_sign =
          xy_displacement_derivative_sign;
    }
    return transits;
  }

  static Time Error(MeasuredTransitsByPlanet const& observations,
                    TransitsByPlanet const& computations,
                    bool const verbose) {
    Exponentiation<Time, 2> sum_error²;
    Time max_error;
    int number_of_transits = 0;
    for (auto const& pair : observations) {
      auto const& name = pair.first;
      auto const& observed_transits = pair.second;
      auto const& computed_transits = computations.at(name);
      for (auto const& observed_transit : observed_transits) {
        auto const next_computed_transit =
            std::lower_bound(computed_transits.begin(),
                             computed_transits.end(),
                             observed_transit.estimated_value);
        Time error;
        if (next_computed_transit == computed_transits.begin()) {
          error = *next_computed_transit - observed_transit.estimated_value;
        } else if (next_computed_transit == computed_transits.end()) {
          error = observed_transit.estimated_value - computed_transits.back();
        } else {
          error = std::min(
              *next_computed_transit - observed_transit.estimated_value,
              observed_transit.estimated_value -
                  *std::prev(next_computed_transit));
        }
        CHECK_LE(0.0 * Second, error);
        LOG_IF(ERROR, verbose)<<name<<": "<<error;
        if (error > max_error) {
          max_error = error;
          LOG_IF(ERROR, verbose)
              << name << ": " 
              << ShortDays(*std::prev(next_computed_transit)) << " "
              << ShortDays(observed_transit.estimated_value) << " "
              << ShortDays(*next_computed_transit) << " " << error;
        }
        sum_error² += Pow<2>(error);
      }
      number_of_transits += observed_transits.size();
    }
    auto const result = Sqrt(sum_error² / number_of_transits);
    LOG_IF(ERROR, verbose)<<"Overall: "<<result<<" "<<number_of_transits;
    return result;
  }

  static std::string SanitizedName(MassiveBody const& body) {
    auto sanitized_name = body.name();
    return sanitized_name.erase(sanitized_name.find_first_of("-"), 1);
  }

  constexpr static char star_name[] = "Trappist-1A";
  SolarSystem<Trappist> const system_;
  not_null<std::unique_ptr<Ephemeris<Trappist>>> ephemeris_;
};

constexpr char TrappistDynamicsTest::star_name[];

TEST_F(TrappistDynamicsTest, MathematicaPeriods) {
  Instant const a_century_later = system_.epoch() + 100 * JulianYear;
  ephemeris_->Prolong(a_century_later);

  auto const& star = system_.massive_body(*ephemeris_, star_name);
  auto const& star_trajectory = ephemeris_->trajectory(star);

  OFStream file(TEMP_DIR / "trappist_periods.generated.wl");
  auto const bodies = ephemeris_->bodies();
  for (auto const& planet : bodies) {
    if (planet != star) {
      auto const& planet_trajectory = ephemeris_->trajectory(planet);
      std::vector<Time> periods;
      for (Instant t = ephemeris_->t_max() - 2000 * Hour;
           t < ephemeris_->t_max();
           t += 1 * Hour) {
        KeplerOrbit<Trappist> const planet_orbit(
            *star,
            *planet,
            planet_trajectory->EvaluateDegreesOfFreedom(t) -
                star_trajectory->EvaluateDegreesOfFreedom(t),
            t);
        periods.push_back(*planet_orbit.elements_at_epoch().period);
      }

      file << mathematica::Assign("period" + SanitizedName(*planet),
                                  periods);
    }
  }
}

TEST_F(TrappistDynamicsTest, MathematicaTransits) {
  Instant const a_century_later = system_.epoch() + 100 * JulianYear;
  ephemeris_->Prolong(a_century_later);

  TransitsByPlanet computations;
  OFStream file(TEMP_DIR / "trappist_transits.generated.wl");

  auto const& star = system_.massive_body(*ephemeris_, star_name);
  auto const bodies = ephemeris_->bodies();
  for (auto const& planet : bodies) {
    if (planet != star) {
      computations[planet->name()] = ComputeTransits(*ephemeris_, star, planet);
      file << mathematica::Assign("transit" + SanitizedName(*planet),
                                  computations[planet->name()]);
    }
  }

  LOG(ERROR) << "max error: "
             << Error(observations, computations, /*verbose=*/true);
  LOG(ERROR) << "χ²: "
             << Transitsχ²(observations, computations, /*verbose=*/true);
}

TEST_F(TrappistDynamicsTest, Optimisation) {
  SolarSystem<Trappist> const system(
      SOLUTION_DIR / "astronomy" / "trappist_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "trappist_initial_state_jd_2457000_000000000.proto.txt");

  bool verbose = false;
  auto planet_names = system.names();
  planet_names.erase(
      std::find(planet_names.begin(), planet_names.end(), star_name));
  std::vector<KeplerianElements<Trappist>> original_elements;
  for (auto const& planet_name : planet_names) {
    original_elements.push_back(SolarSystem<Trappist>::MakeKeplerianElements(
        system.keplerian_initial_state_message(planet_name).elements()));
  }

  double τ = 0.0;
  auto log_pdf_of_system_parameters = 
      [&original_elements, &planet_names, &system, &verbose, &τ](
          SystemParameters const& system_parameters) {
        auto modified_system = system;
        for (int i = 0; i < planet_names.size(); ++i) {
          modified_system.ReplaceElements(
              planet_names[i],
              MakeKeplerianElements(original_elements[i],
                                    system_parameters[i]));
        }

        auto const ephemeris = modified_system.MakeEphemeris(
            /*fitting_tolerance=*/5 * Metre,
            Ephemeris<Trappist>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                                   Position<Trappist>>(),
                /*step=*/0.07 * Day));
        ephemeris->Prolong(modified_system.epoch() + 1000 * Day);

        TransitsByPlanet computations;
        auto const& star = modified_system.massive_body(*ephemeris, star_name);
        auto const bodies = ephemeris->bodies();
        for (auto const& planet : bodies) {
          if (planet != star) {
            computations[planet->name()] =
                ComputeTransits(*ephemeris, star, planet);
          }
        }
        return -Transitsχ²(observations, computations, verbose) / (2.0 * τ);
      };

  Population great_old_ones;
  std::mt19937_64 engine;
  for (int i = 0; i < 50; ++i) {
    great_old_ones.emplace_back();
    SystemParameters& great_old_one = great_old_ones.back();
    std::normal_distribution<> angle_distribution(0.0, 70.0);
    std::normal_distribution<> period_distribution(0.0, 1.0);
    std::normal_distribution<> eccentricity_distribution(0.0, 6.0e-3);
    for (int j = 0; j < great_old_one.size(); ++j) {
      auto perturbed_elements = original_elements[j];
      *perturbed_elements.period += period_distribution(engine) * Second;
      *perturbed_elements.argument_of_periapsis +=
          angle_distribution(engine) * Degree;
      *perturbed_elements.mean_anomaly +=
          angle_distribution(engine) * Degree;
      *perturbed_elements.eccentricity += eccentricity_distribution(engine);
      great_old_one[j] = MakePlanetParameters(perturbed_elements);
    }
  }

  τ = 1000.0;
  RunDEMCMC(great_old_ones,
            /*number_of_generations=*/100,
            /*number_of_generations_between_kicks=*/10,
            /*number_of_burn_in_generations=*/10,
            /*ε=*/0.05,
            log_pdf_of_system_parameters);
}

}  // namespace astronomy
}  // namespace principia
