
#include <memory>

#include "astronomy/frames.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "ksp_plugin/integrators.hpp"
#include "physics/apsides.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"

namespace principia {
namespace astronomy {

using base::not_null;
using geometry::Displacement;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::ComputeApsides;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::MassiveBody;
using physics::SolarSystem;
using quantities::astronomy::JulianYear;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;
using ::testing::Eq;

class SpacecraftTrajectoriesTest : public ::testing::Test {
 protected:
  SpacecraftTrajectoriesTest()
      : solar_system_(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2436145_604166667.proto.txt"),
        ephemeris_(solar_system_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                   Position<ICRS>>(),
                /*step=*/10 * Minute))),
        earth_(solar_system_.massive_body(*ephemeris_, "Earth")),
        jupiter_(solar_system_.massive_body(*ephemeris_, "Jupiter")),
        saturn_(solar_system_.massive_body(*ephemeris_, "Saturn")),
        neptune_(solar_system_.massive_body(*ephemeris_, "Neptune")),
        uranus_(solar_system_.massive_body(*ephemeris_, "Uranus")) {
    google::LogToStderr();
  }

  SolarSystem<ICRS> solar_system_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> ephemeris_;
  not_null<MassiveBody const*> const earth_;
  not_null<MassiveBody const*> const jupiter_;
  not_null<MassiveBody const*> const saturn_;
  not_null<MassiveBody const*> const neptune_;
  not_null<MassiveBody const*> const uranus_;
};

TEST_F(SpacecraftTrajectoriesTest, Voyager2) {
  // An inert Voyager 2.  The flybys are all highly inaccurate, but it should be
  // representative of that kind of trajectory.

  for (int year = 1958; year <= 1977; ++year) {
    Instant const t = ParseUT1(absl::StrCat(year, "-01-01T00:00:00"));
    ephemeris_->Prolong(t);
    LOG(INFO) << "Prolonged to " << year;
  }
  DiscreteTrajectory<ICRS> voyager_2;
  voyager_2.SetDownsampling(ksp_plugin::MaxDenseIntervals,
                            ksp_plugin::DownsamplingTolerance);
  // Data from HORIZONS, identifying TT=TDB.
  voyager_2.Append(
      "1977-08-20T15:33:00"_TT,
      {ICRS::origin +
           Displacement<ICRS>({1.284223722134334E+08 * Kilo(Metre),
                               -7.457444989124168E+07 * Kilo(Metre),
                               -3.234322839044516E+07 * Kilo(Metre)}),
       Velocity<ICRS>({2.151284384548480E+01 * Kilo(Metre) / Second,
                       3.187534896515152E+01 * Kilo(Metre) / Second,
                       1.950130461200638E+01 * Kilo(Metre) / Second})});
  ephemeris_->Prolong(voyager_2.back().time);
  LOG(INFO) << (ephemeris_->trajectory(earth_)->EvaluatePosition(
                    voyager_2.back().time) -
                voyager_2.back().degrees_of_freedom.position())
                       .Norm() /
                   Kilo(Metre)
            << " km from Earth";
  auto const instance = ephemeris_->NewInstance(
      {&voyager_2},
      Ephemeris<ICRS>::NoIntrinsicAccelerations,
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                             Position<ICRS>>(),
          10 * Second));

  for (int year = 1978; year <= 2021; ++year) {
    Instant const last_year = voyager_2.back().time;
    Instant const t = ParseUTC(absl::StrCat(year, "-01-01T00:00:00"));
    ephemeris_->FlowWithFixedStep(t, *instance);
    serialization::DiscreteTrajectory message;
    voyager_2.WriteToMessage(&message, {});
    double const size = 4 * (message.ByteSizeLong() / 3 + 1) / 0x1p10;
    for (auto const body : {jupiter_, saturn_, neptune_, uranus_}) {
      DiscreteTrajectory<ICRS> periapsides;
      DiscreteTrajectory<ICRS> apoapsides;
      ComputeApsides(*ephemeris_->trajectory(body),
                     voyager_2.LowerBound(last_year),
                     voyager_2.end(),
                     /*max_points=*/std::numeric_limits<int>::max(),
                     apoapsides,
                     periapsides);
      if (!periapsides.Empty()) {
        LOG(INFO) << (periapsides.back().degrees_of_freedom.position() -
                      ephemeris_->trajectory(body)->EvaluatePosition(
                          periapsides.back().time))
                             .Norm() /
                         Kilo(Metre)
                  << " km from " << body->name();
      }
    }
    LOG(INFO) << "Happy new year " << year << "; trajectory has "
              << voyager_2.Size() << " points, " << size << " kiB, "
              << size / (t - voyager_2.front().time) * JulianYear << " kiB/a";
  }

  // From the IMCCE calendar for 2021.
  constexpr Instant gröbner_release_date = "2021-06-10T10:52:39"_UTC;
  ephemeris_->FlowWithFixedStep(gröbner_release_date, *instance);
}

}  // namespace astronomy
}  // namespace principia