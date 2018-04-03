
#include <iostream>
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
    if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " cbrt (atlas|egg|kahan|microsoft)";
    }
    std::string method = argv[2];
    std::cout << "{";
    std::uint64_t Y;
    std::cin >> Y;
    for (;;) {
      double y = principia::to_double(Y);
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
      std::cout << principia::to_integer(x);
      std::cin >> Y;
      if (Y == 0) {
        std::cout << "}\n";
        break;
      } else {
        std::cout << ",";
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
