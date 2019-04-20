
#include <limits>
#include <random>

#include "absl/strings/str_cat.h"
#include "astronomy/standard_product_3.hpp"
#include "base/bundle.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"

namespace principia {
namespace astronomy {

using base::Bundle;
using base::not_null;
using base::Status;
using geometry::AngleBetween;
using geometry::Bivector;
using geometry::Displacement;
using geometry::Instant;
using geometry::Normalize;
using geometry::OrientedAngleBetween;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;
using geometry::Velocity;
using geometry::Wedge;
using integrators::EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Fine1987RKNG34;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodySurfaceDynamicFrame;
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::Acceleration;
using quantities::Angle;
using quantities::ArcSin;
using quantities::Cos;
using quantities::Length;
using quantities::Pow;
using quantities::Sin;
using quantities::Sqrt;
using quantities::Square;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::testing::TestWithParam;
using ::testing::ValuesIn;

namespace {

struct StandardProduct3Args {
  std::filesystem::path filename;
  StandardProduct3::Dialect dialect = StandardProduct3::Dialect::Standard;
};

std::ostream& operator<<(std::ostream& out, StandardProduct3Args const& args) {
  return out << args.filename << " interpreted as " << args.dialect;
}

namespace deмcmc {

// TODO(egg): We should factor the two implementations here and in the
// TRAPPIST-1 optimization, and move that to a dedicated file.

template<typename Parameters>
std::vector<double> EvaluatePopulation(
    std::vector<Parameters> const& population,
    std::function<double(Parameters const&, std::string&)> const& compute_log_pdf,
    std::vector<std::string>& info) {
  std::vector<double> log_pdf(population.size());
  info.resize(population.size());
  Bundle bundle;
  for (int i = 0; i < population.size(); ++i) {
    auto const& parameters = population[i];
    bundle.Add([&compute_log_pdf, i, &log_pdf, &parameters, &info]() {
      log_pdf[i] = compute_log_pdf(parameters, info[i]);
      return Status::OK;
    });
  }
  bundle.Join();
  return log_pdf;
}

template<typename Parameters>
std::vector<Parameters> GenerateTrialStates(
    std::vector<Parameters> const& population,
    double const γ,
    double const ε,
    std::mt19937_64& engine) {
  std::vector<Parameters> trial(population.size());
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

template<typename Parameters>
Parameters Run(std::vector<Parameters>& population,
               int const number_of_generations,
               int const number_of_generations_between_kicks,
               int const number_of_burn_in_generations,
               double const ε,
               std::function<double(Parameters const&, std::string&)> const&
                   compute_log_pdf) {
  CHECK_GE(number_of_generations, 1);
  std::mt19937_64 engine;
  std::uniform_real_distribution<> distribution(0.0, 1.0);

  static double best_log_pdf = -std::numeric_limits<double>::infinity();
  std::string best_info;
  Parameters best_parameters;

  std::vector<std::string> infos;
  auto log_pdf = EvaluatePopulation(population, compute_log_pdf, infos);

  // Loop over generations.
  for (int generation = 0; generation < number_of_generations; ++generation) {
    int accepted = 0;

    // Every 10th generation try full-size steps.
    double const γ = generation < number_of_burn_in_generations
                         ? 0.01
                         : generation % number_of_generations_between_kicks == 0
                               ? 1.0
                               : 2.38 / Sqrt(2 * 11);

    // Evaluate model for each set of trial parameters.
    auto const trial = GenerateTrialStates(population, γ, ε, engine);
    std::vector<std::string> trial_infos;
    auto const log_pdf_trial =
        EvaluatePopulation(trial, compute_log_pdf, trial_infos);

    // For each member of population.
    for (int i = 0; i < population.size(); ++i) {
      double const log_pdf_ratio = log_pdf_trial[i] - log_pdf[i];
      if (log_pdf_ratio > 0.0 ||
          log_pdf_ratio > std::log(distribution(engine))) {
        population[i] = trial[i];
        log_pdf[i] = log_pdf_trial[i];
        infos[i] = trial_infos[i];
        ++accepted;
      }
    }

    // Traces.
    int const max_index =
        std::max_element(log_pdf.begin(), log_pdf.end()) - log_pdf.begin();
    if (best_log_pdf < log_pdf[max_index]) {
      best_parameters = population[max_index];
      best_log_pdf = log_pdf[max_index];
      best_info = infos[max_index];
    }
    LOG(ERROR) << "Generation " << generation << "; Acceptance: " << accepted
               << " / " << population.size();
    LOG(ERROR) << "Max  : " << infos[max_index];
    if (best_info != infos[max_index]) {
      LOG(ERROR) << "Best : " << best_info;
    }
  }
  return best_parameters;
}

}  // namespace deмcmc

}  // namespace

class SolarRadiationPressureTest : public TestWithParam<StandardProduct3Args> {
 public:
  static void SetUpTestCase() {
    google::LogToStderr();
    if (static_ephemeris_ == nullptr) {
      solar_system_2010_.LimitOblatenessToDegree("Earth", 2);
      static_ephemeris_ = solar_system_2010_.MakeEphemeris(
          /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                   /*geopotential_tolerance=*/0x1p-24},
          Ephemeris<ICRS>::FixedStepParameters(
              SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                 Position<ICRS>>(),
              /*step=*/10 * Minute));
      for (int year = 2010; year <= 2019; ++year) {
        std::string date = absl::StrCat(year, "-06-01T00:00:00");
        LOG(INFO) << "Prolonging to " << date << " (UTC)";
        static_ephemeris_->Prolong(ParseUTC(date));
      }
    }
  }

 protected:
  SolarRadiationPressureTest()
      : sp3_(GetParam().filename, GetParam().dialect),
        ephemeris_(*static_ephemeris_),
        earth_(solar_system_2010_.rotating_body(ephemeris_, "Earth")),
        sun_(solar_system_2010_.rotating_body(ephemeris_, "Sun")),
        earth_trajectory_(*ephemeris_.trajectory(earth_)),
        sun_trajectory_(*ephemeris_.trajectory(sun_)),
        itrs_(&ephemeris_, earth_) {}

  double Computeχ²(StandardProduct3::SatelliteIdentifier const& satellite,
                   physics::RelativeDegreesOfFreedom<ITRS> const& initial_velocity,
                   Ephemeris<ICRS>::GeneralizedIntrinsicAcceleration const&
                       solar_radiation_pressure,
                   std::string& info) const {
    auto const initial_it = sp3_.orbit(satellite).front()->Begin();
    DiscreteTrajectory<ICRS> integrated;
    integrated.Append(initial_it.time(),
                      itrs_.FromThisFrameAtTime(initial_it.time())(
                          initial_it.degrees_of_freedom() + initial_velocity));
    double χ² = 0;
    Square<Length> max_residual;
    int n = 0;
    for (not_null<DiscreteTrajectory<ITRS> const*> const arc :
         sp3_.orbit(satellite)) {
      for (auto it = arc->Begin(); it != arc->End(); ++it, ++n) {
        if (it.time() < initial_it.time()) {
          continue;
        }
        Ephemeris<ICRS>::GeneralizedAdaptiveStepParameters parameters(
            EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
                Fine1987RKNG34,
                Position<ICRS>>(),
            /*max_steps=*/std::numeric_limits<int64_t>::max(),
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Micro(Metre) / Second);
        ephemeris_.FlowWithAdaptiveStep(
            &integrated,
            solar_radiation_pressure,
            it.time(),
            parameters,
            /*max_ephemeris_steps=*/std::numeric_limits<int64_t>::max(),
            /*last_point_only=*/true);
        Square<Length> const residual =
            (itrs_.ToThisFrameAtTime(it.time())(
                 integrated.last().degrees_of_freedom()).position() -
             it.degrees_of_freedom().position()).Norm²();
        // TODO(egg): use an actual uncertainty here.
        χ² += residual / Pow<2>(Pow<5>(1.25) * Milli(Metre));
        max_residual = std::max(max_residual, residual);
      }
    }
    info = (std::stringstream() << u8"χ² = " << χ²  << u8" χ²/ν = " <<
            χ² / (3 * (integrated.Size() - 1) - 11) << " max error = "
            << Sqrt(max_residual)).str();
    return χ²;
  }

 private:
  static SolarSystem<ICRS> solar_system_2010_;
  static std::unique_ptr<Ephemeris<ICRS>> static_ephemeris_;

 protected:
  StandardProduct3 const sp3_;
  Ephemeris<ICRS>& ephemeris_;
  not_null<RotatingBody<ICRS> const*> const earth_;
  not_null<RotatingBody<ICRS> const*> const sun_;
  ContinuousTrajectory<ICRS> const& earth_trajectory_;
  ContinuousTrajectory<ICRS> const& sun_trajectory_;
  BodySurfaceDynamicFrame<ICRS, ITRS> const itrs_;
};

SolarSystem<ICRS> SolarRadiationPressureTest::solar_system_2010_(
    SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
    SOLUTION_DIR / "astronomy" /
        "sol_initial_state_jd_2455200_500000000.proto.txt");
std::unique_ptr<Ephemeris<ICRS>> SolarRadiationPressureTest::static_ephemeris_;

INSTANTIATE_TEST_CASE_P(NGA,
                        SolarRadiationPressureTest,
                        ValuesIn(std::vector<StandardProduct3Args>{
                            {SOLUTION_DIR / "astronomy" / "standard_product_3" /
                                 "COD0MGXFIN_20181260000_01D_05M_ORB.SP3",
                             StandardProduct3::Dialect::ChineseMGEX},
                            {SOLUTION_DIR / "astronomy" / "standard_product_3" /
                             "COD0MGXFIN_20181260000_01D_05M_ORB.SP3"},
                        }));

struct DYB;
struct ReducedECOMParameters {
  physics::RelativeDegreesOfFreedom<ITRS> initial_velocity = {
      Displacement<ITRS>(),
      Velocity<ITRS>()};
  Vector<Acceleration, DYB> order_0_acceleration;
  Acceleration Bc;
  Acceleration Bs;
};

std::ostream& operator<<(std::ostream& out, ReducedECOMParameters const& parameters) {
  return out << parameters.order_0_acceleration << ", (" << parameters.Bc
             << ", " << parameters.Bs << "), " << parameters.initial_velocity;
}

ReducedECOMParameters operator+(ReducedECOMParameters const& left,
                                ReducedECOMParameters const& right) {
  return {left.initial_velocity + right.initial_velocity,
          left.order_0_acceleration + right.order_0_acceleration,
          left.Bc + right.Bc,
          left.Bs + right.Bs};
}

ReducedECOMParameters operator-(ReducedECOMParameters const& left,
                                ReducedECOMParameters const& right) {
  return {left.initial_velocity - right.initial_velocity,
          left.order_0_acceleration - right.order_0_acceleration,
          left.Bc - right.Bc,
          left.Bs - right.Bs};
}

ReducedECOMParameters operator*(double const left,
                                ReducedECOMParameters const& right) {
  return {left * right.initial_velocity,
          left * right.order_0_acceleration,
          left * right.Bc,
          left * right.Bs};
}

TEST_P(SolarRadiationPressureTest, ReducedECOM) {/*
  LOG(ERROR) << "===============================================";
  Instant const t = "2010-06-01T12:00:00"_TT;
  LOG(ERROR) << sun_trajectory_.EvaluateDegreesOfFreedom(t);
  LOG(ERROR) << sun_trajectory_.EvaluateDegreesOfFreedom(t) -
                    earth_trajectory_.EvaluateDegreesOfFreedom(t);
  LOG(ERROR) << itrs_.ToThisFrameAtTime(t)(
      sun_trajectory_.EvaluateDegreesOfFreedom(t));
  LOG(ERROR) << (itrs_
                     .ToThisFrameAtTime(t)(
                         sun_trajectory_.EvaluateDegreesOfFreedom(t))
                     .position() -
                 ITRS::origin)
                    .coordinates()
                    .ToSpherical();
  LOG(ERROR) << itrs_
                    .ToThisFrameAtTime(t)(
                        sun_trajectory_.EvaluateDegreesOfFreedom(t))
                    .velocity()
                    .coordinates()
                    .ToSpherical();
  LOG(FATAL) << "===============================================";*/
  for (auto const& satellite : sp3_.satellites()) {
    ReducedECOMParameters candidate;
    std::mt19937_64 engine(1729);
    auto solar_radiation_pressure =
        [this](ReducedECOMParameters const parameters,
               Instant const& t,
               DegreesOfFreedom<ICRS> const& satellite_dof)
        -> Vector<Acceleration, ICRS> {
      DegreesOfFreedom<ICRS> const earth_dof =
          earth_trajectory_.EvaluateDegreesOfFreedom(t);
      Position<ICRS> const sun_position = sun_trajectory_.EvaluatePosition(t);
      Displacement<ICRS> const satellite_to_sun =
          sun_position - satellite_dof.position();
      // Notation from Arnold et al. (2015), CODE’s new solar radiation pressure
      // model for GNSS orbit determination.
      Displacement<ICRS> const r = satellite_dof.position() - earth_dof.position();
      Vector<double, ICRS> const eD = Normalize(satellite_to_sun);
      Vector<double, ICRS> const er = Normalize(r);
      Bivector<double, ICRS> const eY = -Normalize(Wedge(er, eD));
      Vector<double, ICRS> const eB = eD * eY;
      Rotation<DYB, ICRS> const to_icrs(eD, eY, eB);

      Velocity<ICRS> const v = satellite_dof.velocity() - earth_dof.velocity();
      Bivector<double, ICRS> const orbit_normal = Normalize(Wedge(r, v));
      Vector<double, ICRS> const north({0, 0, 1});
      Vector<double, ICRS> const ascending_node = north * orbit_normal;

      Angle const u = OrientedAngleBetween(ascending_node, r, orbit_normal);
      // Shadow of a point sun; alternatively, we could rigorously compute the
      // umbra and penumbra.
      if (AngleBetween(satellite_to_sun, -r) <
          ArcSin(earth_->mean_radius() / r.Norm())) {
        return Vector<Acceleration, ICRS>();
      }
      return to_icrs(parameters.order_0_acceleration +
                     Vector<Acceleration, DYB>(
                         {0 * Metre / Pow<2>(Second),
                          0 * Metre / Pow<2>(Second),
                          parameters.Bc * Cos(u) + parameters.Bs * Sin(u)}));
    };
    std::string info;
    Computeχ²(satellite,
              {Displacement<ITRS>{}, Velocity<ITRS>{}},
              std::bind(solar_radiation_pressure,
                        ReducedECOMParameters{}, _1, _2),
              info);
    std::cout << satellite << " no SRP: " << info << "\n";
    std::vector<ReducedECOMParameters> population(
        120,
        ReducedECOMParameters{
            {Displacement<ITRS>{}, Velocity<ITRS>{}},
            Vector<Acceleration, DYB>({-1e-7 * Metre / Pow<2>(Second),
                                       0 * Metre / Pow<2>(Second),
                                       0 * Metre / Pow<2>(Second)}),
            0 * Metre / Pow<2>(Second),
            0 * Metre / Pow<2>(Second)});
    for (auto& parameters : population) {
      parameters.initial_velocity +=
          {Displacement<ITRS>(
               {std::normal_distribution()(engine) * (1e-3 * Metre),
                std::normal_distribution()(engine) * (1e-3 * Metre),
                std::normal_distribution()(engine) * (1e-3 * Metre)}),
           Velocity<ITRS>(
               {std::normal_distribution()(engine) * (1e-4 * Metre / Second),
                std::normal_distribution()(engine) * (1e-4 * Metre / Second),
                std::normal_distribution()(engine) * (1e-4 * Metre / Second)})};
      parameters.order_0_acceleration += Vector<Acceleration, DYB>(
          {std::normal_distribution()(engine) * (1e-8 * Metre / Pow<2>(Second)),
           std::normal_distribution()(engine) * (1e-8 * Metre / Pow<2>(Second)),
           std::normal_distribution()(engine) *
               (1e-8 * Metre / Pow<2>(Second))});
      parameters.Bc +=
          std::normal_distribution()(engine) * (1e-8 * Metre / Pow<2>(Second));
      parameters.Bs +=
          std::normal_distribution()(engine) * (1e-8 * Metre / Pow<2>(Second));
    }
    deмcmc::Run<ReducedECOMParameters>(
        population,
        1'147,
        30,
        10,
        0.05,
        [this, satellite, &solar_radiation_pressure](
            ReducedECOMParameters const& parameters,
            std::string& info) {
          double const χ² =
              Computeχ²(
                  satellite,
                  parameters.initial_velocity,
                  std::bind(solar_radiation_pressure, parameters, _1, _2),
                  info);
          info += (std::stringstream() << ": " << parameters).str();
          return -χ² / 2;
        });
        return;
  }
}

}  // namespace astronomy
}  // namespace principia
