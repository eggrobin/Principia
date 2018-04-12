
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
    std::int64_t iterations = std::atoll(argv[3]);
    std::mt19937_64 mersenne;
    std::uint64_t const binary64_1 = principia::to_integer(1);
    std::uint64_t const binary64_8 = principia::to_integer(8);
    int incorrect_roundings = 0;
    int unfaithful_roundings = 0;
    double max_ulps = 0;
    double min_positive_correct_ulps =
        std::numeric_limits<double>::infinity();
    std::uint64_t Y_just_below_exact;
    double max_negative_correct_ulps =
        -std::numeric_limits<double>::infinity();
    std::uint64_t Y_just_above_exact;
    double min_positive_incorrect_ulps =
        std::numeric_limits<double>::infinity();
    std::uint64_t Y_just_below_tie;
    double max_negative_incorrect_ulps =
        -std::numeric_limits<double>::infinity();
    std::uint64_t Y_just_above_tie;
    std::uint64_t Y_worst;
    for (std::int64_t i = 1; i <= iterations; ++i) {
      std::uint64_t const Y =
          mersenne() % (binary64_8 - binary64_1) + binary64_1;
      double const y = principia::to_double(Y);
      principia::slow_correct::RoundedReal const x_correct =
          principia::slow_correct::cube_root(y);
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
      } else if (method == "sun") {
        x = principia::sun::cbrt(y);
      }
      std::int64_t const ulps_from_correct =
          principia::to_integer(x) -
          principia::to_integer(x_correct.nearest_rounding);
      double const ulps_from_exact = x_correct.nearest_ulps + ulps_from_correct;
      double const abs_ulps_from_exact = std::abs(ulps_from_exact);
      if (abs_ulps_from_exact > max_ulps) {
        max_ulps = abs_ulps_from_exact;
        Y_worst = Y;
      }
      if (x != x_correct.nearest_rounding) {
        ++incorrect_roundings;
        CHECK_GT(abs_ulps_from_exact, 0.5);
        if (x != x_correct.furthest_rounding) {
          ++unfaithful_roundings;
          CHECK_GT(abs_ulps_from_exact, 1);
        }
      } else {
        CHECK_LE(abs_ulps_from_exact, 0.5);
      }

      std::int64_t const ulps_incorrect_from_correct =
          principia::to_integer(x_correct.furthest_rounding) -
          principia::to_integer(x_correct.nearest_rounding);
      double const ulps_incorrect_from_exact =
          x_correct.nearest_ulps + ulps_incorrect_from_correct;
      if (ulps_incorrect_from_exact > 0) {
        if (ulps_incorrect_from_exact < min_positive_incorrect_ulps) {
          min_positive_incorrect_ulps = ulps_incorrect_from_exact;
          Y_just_below_tie = Y;
        }
      } else {
        if (ulps_incorrect_from_exact > max_negative_incorrect_ulps) {
          max_negative_incorrect_ulps = ulps_incorrect_from_exact;
          Y_just_above_tie = Y;
        }
      }
      if (x_correct.nearest_ulps != 0) {
        if (x_correct.nearest_ulps > 0) {
          if (x_correct.nearest_ulps < min_positive_correct_ulps) {
            min_positive_correct_ulps = x_correct.nearest_ulps;
            Y_just_below_exact = Y;
          }
        } else {
          if (x_correct.nearest_ulps > max_negative_correct_ulps) {
            max_negative_correct_ulps = x_correct.nearest_ulps;
            Y_just_above_exact = Y;
          }
        }
      }

      if (i % 1'000'000 == 0) {
        std::cout << "Method " << method << ".\n"
                  << "Tested " << i << " values in [1, 8[.\n"
                  << "incorrect roundings  : " << incorrect_roundings << " ("
                  << 100.0 * incorrect_roundings / i << " %)\n"
                  << "unfaithful roundings : " << unfaithful_roundings << " ("
                  << 100.0 * unfaithful_roundings / i << " %)\n"
                  << "maximal error        : " << max_ulps << " ULPs, for "
                  << std::hex << std::uppercase << "16^^" << Y_worst << std::dec
                  << "\n\n"
                  << "Trickiest sample values:\n"
                  << "closest to ties      : 1/2 + "
                  << min_positive_incorrect_ulps - 0.5 << " ULPs, for "
                  << std::hex << std::uppercase << "16^^" << Y_just_below_tie
                  << std::dec << "\n"
                  << "                      -1/2 - "
                  << -0.5 - max_negative_incorrect_ulps << " ULPs, for "
                  << std::hex << std::uppercase << "16^^" << Y_just_above_tie
                  << std::dec << "\n"
                  << "closest to exact     :       "
                  << min_positive_correct_ulps << " ULPs, for " << std::hex
                  << std::uppercase << "16^^" << Y_just_below_exact << std::dec
                  << "\n"
                  << "                            " << max_negative_correct_ulps
                  << " ULPs, for " << std::hex << std::uppercase << "16^^"
                  << Y_just_above_exact << std::dec << "\n";
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
