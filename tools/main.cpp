
#include <iostream>
#include <random>
#include <string>

#include "astronomy/epoch.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "quantities/parser.hpp"
#include "tools/generate_configuration.hpp"
#include "tools/generate_profiles.hpp"

#define NO_BENCHMARK 1
#include "benchmarks/cbrt.cpp"

int main(int argc, char const* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::LogToStderr();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " command [arguments...]\n";
    return 1;
  }
  std::string command = argv[1];
  if (command == "cbrt") {
    if (argc != 4) {
      std::cerr << "Usage: " << argv[0]
                << " cbrt (atlas|egg|kahan|microsoft) iterations";
    }
    std::string method = argv[2];
    int iterations = std::atoi(argv[3]);
    std::mt19937_64 mersenne;
    std::uint64_t const binary64_1 = principia::to_integer(1);
    std::uint64_t const binary64_8 = principia::to_integer(8);
    int incorrect_roundings = 0;
    int unfaithful_roundings = 0;
    for (int i = 1; i <= iterations; ++i) {
      std::uint64_t const Y =
          mersenne() % (binary64_8 - binary64_1) + binary64_1;
      double const y = principia::to_double(Y);
      double x_incorrect;
      double const x_correct = principia::slow_correct::cbrt(y, &x_incorrect);
      double x;
      if (method == "atlas") {
        x = principia::atlas::cbrt(y);
      } else if (method == "householder_order_10_estrin") {
        x = principia::householder_order_10_estrin::cbrt(y);
      } else if (method == "egg") {
        x = principia::egg::cbrt(y);
      } else if (method == "kahan") {
        x = principia::kahan::cbrt(y);
      } else if (method == "microsoft") {
        x = std::cbrt(y);
      }
      if (x != x_correct) {
        ++incorrect_roundings;
        if (x != x_incorrect) {
          ++unfaithful_roundings;
        }
      }
      if (i % 1'000'000 == 0) {
        std::cout << "Tested " << i << " values in [1, 8[.\n"
                  << "incorrect roundings  : " << incorrect_roundings << "("
                  << 100.0 * incorrect_roundings / i << " %)\n"
                  << "unfaithful roundings : " << unfaithful_roundings << "("
                  << 100.0 * unfaithful_roundings / i << " %)\n";
      }
    }
  } else if (command == "generate_configuration") {
    if (argc != 6) {
      // tools.exe generate_configuration \
      //     JD2433647.5 \
      //     sol_gravity_model \
      //     sol_initial_state_jd_2433282_500000000 \
      //     sol_numerics_blueprint
      std::cerr << "Usage: " << argv[0] << " " << argv[1] << " "
                << "game_epoch_jd "
                << "gravity_model_stem "
                << "initial_state_stem "
                << "numerics_blueprint_stem\n";
      return 2;
    }
    std::string const game_epoch = argv[2];
    std::string const gravity_model_stem = argv[3];
    std::string const initial_state_stem = argv[4];
    std::string const numerics_blueprint_stem = argv[5];
    principia::tools::GenerateConfiguration(game_epoch,
                                            gravity_model_stem,
                                            initial_state_stem,
                                            numerics_blueprint_stem);
    return 0;
  } else if (command == "generate_profiles") {
    if (argc != 2) {
      // tools.exe generate_profiles
      std::cerr << "Usage: " << argv[0] << " " << argv[1] << "\n";
      return 3;
    }
    principia::tools::GenerateProfiles();
    return 0;
  } else {
    std::cerr << "Usage: " << argv[0]
              << " generate_configuration|generate_profiles\n";
    return 4;
  }
}
