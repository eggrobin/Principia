
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "absl/strings/str_replace.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_cat.h"
#include "astronomy/epoch.hpp"
#include "benchmarks/cbrt.hpp"
#include "geometry/named_quantities.hpp"
#include "glog/logging.h"
#include "quantities/parser.hpp"
#include "tools/generate_configuration.hpp"
#include "tools/generate_kopernicus.hpp"
#include "tools/generate_profiles.hpp"

namespace {
std::string const MeasuredProportion(double const successes,
                                     double const samples) {
  if (samples == 0) {
    return "NaN";
  } else if (successes == 0) {
    return "0";
  } else if (successes == samples) {
    return "1";
  }
  double const estimated_proportion = successes / samples;
  double const u = principia::quantities::Sqrt(
      estimated_proportion * (1 - estimated_proportion) / (samples - 1));
  std::int64_t const floor_log10_u = std::floor(std::log10(u));
  std::int64_t digits_shown =
      std::floor(std::log10(estimated_proportion)) - floor_log10_u + 2;
  std::int64_t uncertainty_parenthetical =
      std::ceil(u * std::pow(10, 1 - floor_log10_u));
  return absl::StrReplaceAll(
      absl::StrFormat(absl::StrCat("%1.", digits_shown - 1, "e").c_str(),
                      estimated_proportion),
      {{"e", absl::StrFormat(u8"(%02d)×10^", uncertainty_parenthetical)}});
}
}  // namespace

int __cdecl main(int argc, char const* argv[]) {
  google::SetLogFilenameExtension(".log");
  google::InitGoogleLogging(argv[0]);
  google::LogToStderr();
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " command [arguments...]\n";
    return 1;
  }
  std::string command = argv[1];
  if (command == "cbrt" || command == "cbrtfast") {
    const bool fast = command == "cbrtfast";
    if (argc != 4) {
      std::cerr << "Usage: " << argv[0] << " cbrt[fast] methods iterations";
    }
    std::string const method_arg = argv[2];
    std::int64_t iterations = std::atoll(argv[3]);
    std::mt19937_64 mersenne;
    std::uint64_t const binary64_1 = principia::numerics::to_integer(1);
    std::uint64_t const binary64_8 = principia::numerics::to_integer(8);
    struct MethodProperties {
      int incorrect_roundings = 0;
      int unfaithful_roundings = 0;
      double max_ulps = 0;
      std::uint64_t Y_worst = 0;
      int detected_possible_misroundings = 0;
      int undetected_misroundings = 0;
    };
    std::map<std::string, MethodProperties> properties;
    double min_positive_correct_ulps = std::numeric_limits<double>::infinity();
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
    std::vector<std::string> methods;
    if (method_arg == "all") {
      for (auto const& pair :
           principia::numerics::CubeRootRegistry::Instance().methods()) {
        if (fast && pair.first == "fast_correct") {
          continue;
        }
        methods.emplace_back(pair.first);
      }
    } else {
      for (int begin = 0, end = 0; end != std::string::npos; begin = end + 1) {
        end = method_arg.find(',', begin);
        methods.emplace_back(method_arg.substr(begin, end - begin));
      }
    }
    std::cout << "Methods:\n";
    for (auto const& method : methods) {
      std::cout << "  " << method << "\n";
    }
    for (std::int64_t i = 1; i <= iterations; ++i) {
      std::uint64_t const Y =
          mersenne() % (binary64_8 - binary64_1) + binary64_1;
      double const y = principia::numerics::to_double(Y);
      double x_nearest;
      principia::numerics::RoundedReal x_correct;
      if (fast) {
        x_nearest = principia::numerics::fast_correct::cbrt(y);
      } else {
        x_correct = principia::numerics::correct_cube_root(y);
        x_nearest = x_correct.nearest_rounding;
      }
      for (auto const& method : methods) {
        principia::numerics::possible_misrounding = false;
        double const x =
            principia::numerics::CubeRootRegistry::Instance().methods().at(
                method)(y);
        if (principia::numerics::possible_misrounding) {
          ++properties[method].detected_possible_misroundings;
        }
        std::int64_t const ulps_from_correct =
            principia::numerics::to_integer(x) -
            principia::numerics::to_integer(x_nearest);
        double abs_ulps_from_exact;
        if (!fast) {
          double const ulps_from_exact =
              x_correct.nearest_ulps + ulps_from_correct;
          abs_ulps_from_exact = std::abs(ulps_from_exact);
        }
        if (abs_ulps_from_exact > properties[method].max_ulps) {
          properties[method].max_ulps = abs_ulps_from_exact;
          properties[method].Y_worst = Y;
        }
        if (x != x_nearest) {
          ++properties[method].incorrect_roundings;
          if (!principia::numerics::possible_misrounding) {
            ++properties[method].undetected_misroundings;
          }
          if (!fast) {
            CHECK_GT(abs_ulps_from_exact, 0.5);
            if (x != x_correct.furthest_rounding) {
              ++properties[method].unfaithful_roundings;
              CHECK_GT(abs_ulps_from_exact, 1);
            }
          }
        } else if (!fast) {
          CHECK_LE(abs_ulps_from_exact, 0.5);
        }
      }

      if (!fast) {
        std::int64_t const ulps_incorrect_from_correct =
            principia::numerics::to_integer(x_correct.furthest_rounding) -
            principia::numerics::to_integer(x_nearest);
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
      }

      if (i % 1'000'000 == 0) {
        std::cout << "Tested " << i << " values in [1, 8[.\n\n";
        if (!fast) {
          std::cout << "Trickiest sample values:\n"
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
                    << std::uppercase << "16^^" << Y_just_below_exact
                    << std::dec << "\n"
                    << "                            "
                    << max_negative_correct_ulps << " ULPs, for " << std::hex
                    << std::uppercase << "16^^" << Y_just_above_exact
                    << std::dec << "\n\n";
        }
        for (auto const& method : methods) {
          std::cout << "Method " << method << "\n"
                    << "incorrect roundings             : "
                    << properties[method].incorrect_roundings << ", "
                    << MeasuredProportion(
                           properties[method].incorrect_roundings, i)
                    << "\n";
          if (!fast) {
            std::cout << "unfaithful roundings            : "
                      << properties[method].unfaithful_roundings << ", "
                      << MeasuredProportion(
                             properties[method].unfaithful_roundings, i)
                      << "\n";
          }
          if (properties[method].detected_possible_misroundings != 0) {
            std::cout << "detected possible misroundings  : "
                      << properties[method].detected_possible_misroundings
                      << ", "
                      << MeasuredProportion(
                             properties[method].detected_possible_misroundings,
                             i)
                      << "\n"
                      << "undetected misroundings         : "
                      << properties[method].undetected_misroundings << ", "
                      << MeasuredProportion(
                             properties[method].undetected_misroundings,
                             properties[method].incorrect_roundings)
                      << " of misroundings\n";
          }
          if (!fast) {
            std::cout << "maximal error                   : "
                      << properties[method].max_ulps << " ULPs, for "
                      << std::hex << std::uppercase << "16^^"
                      << properties[method].Y_worst << std::dec << "\n";
          }
          std::cout << "\n";
        }
      }
    }
  } else if (command == "generate_configuration") {
    if (argc != 7) {
      // tools.exe generate_configuration \
      //     JD2433647.5 \
      //     sol_gravity_model \
      //     sol_initial_state_jd_2433282_500000000 \
      //     sol_numerics_blueprint \
      //     RealSolarSystem
      // tools.exe generate_configuration \
      //     JD2457000.0 \
      //     trappist_gravity_model \
      //     trappist_initial_state_jd_2457000_000000000 \
      //     trappist_numerics_blueprint \
      //     aSLIPPIST-1
      std::cerr << "Usage: " << argv[0] << " " << argv[1] << " "
                << "game_epoch_jd "
                << "gravity_model_stem "
                << "initial_state_stem "
                << "numerics_blueprint_stem "
                << "needs\n";
      return 2;
    }
    std::string const game_epoch = argv[2];
    std::string const gravity_model_stem = argv[3];
    std::string const initial_state_stem = argv[4];
    std::string const numerics_blueprint_stem = argv[5];
    std::string const needs = argv[6];
    principia::tools::GenerateConfiguration(game_epoch,
                                            gravity_model_stem,
                                            initial_state_stem,
                                            numerics_blueprint_stem,
                                            needs);
    return 0;
  } else if (command == "generate_kopernicus") {
    if (argc != 4) {
      // tools.exe generate_kopernicus \
      //     trappist_gravity_model \
      //     trappist_initial_state_jd_2457000_000000000
      std::cerr << "Usage: " << argv[0] << " " << argv[1] << " "
                << "gravity_model_stem "
                << "initial_state_stem\n";
      return 5;
    }
    std::string const gravity_model_stem = argv[2];
    std::string const initial_state_stem = argv[3];
    principia::tools::GenerateKopernicusForSlippist1(gravity_model_stem,
                                                     initial_state_stem);
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
