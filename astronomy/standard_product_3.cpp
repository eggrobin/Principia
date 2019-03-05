﻿#include "astronomy/standard_product_3.hpp"

#include <fstream>
#include <optional>

#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "astronomy/time_scales.hpp"
#include "base/map_util.hpp"
#include "glog/logging.h"

namespace principia {
namespace astronomy {
namespace internal_standard_product_3 {

using base::FindOrDie;
using geometry::Displacement;
using quantities::NaN;
using quantities::Speed;
using quantities::si::Deci;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Second;

StandardProduct3::StandardProduct3(
    std::filesystem::path const& filename,
    StandardProduct3::Dialect const dialect) {
  std::ifstream file(filename);
  CHECK(file.good()) << filename;
  std::optional<std::string> line;

  std::string location;
  int line_number = 0;
  auto const read_line = [&file, &filename, &line, &line_number, &location]() {
    if (!line.has_value()) {
      line.emplace();
    }
    std::getline(file, *line);
    if (file.fail()) {
      CHECK(file.eof()) << "Non-EOF failure after " << location;
      line.reset();
      location = absl::StrCat(filename.string(), " at end of file");
    } else {
      ++line_number;
      location =
          absl::StrCat(filename.string(), " line ", line_number, ": ", *line);
    }
  };

  // The specification uses 1-based column indices, and column ranges with
  // bounds included.
  auto const column = [&line, &location, line_number](int const index) {
    CHECK(line.has_value()) << location;
    CHECK_LT(index - 1, line->size()) << location;
    return (*line)[index - 1];
  };
  auto const columns = [&line, &location](int const first, int const last) {
    CHECK(line.has_value()) << location;
    CHECK_LT(last - 1, line->size()) << location;
    CHECK_LE(first, last) << location;
    return std::string_view(&(*line)[first - 1], last - first + 1);
  };
  auto const float_columns = [&columns, &location](int const first,
                                                   int const last) {
    double result;
    CHECK(absl::SimpleAtod(columns(first, last), &result))
        << location << " columns " << first << "-" << last;
    return result;
  };
  auto const integer_columns = [&columns, &location](int const first,
                                                     int const last) {
    int result;
    CHECK(absl::SimpleAtoi(columns(first, last), &result))
        << location << " columns " << first << "-" << last;
    return result;
  };

  int number_of_epochs;
  int number_of_satellites;

  // Header: # record.
  read_line();
  CHECK_EQ(column(1), '#') << location;
  CHECK_GE(Version{column(2)}, Version::A) << location;
  CHECK_LE(Version{column(2)}, Version::D) << location;
  version_ = Version{column(2)};
  CHECK(column(3) == 'P' || column(3) == 'V') << location;
  has_velocities_ = column(3) == 'V';
  number_of_epochs = integer_columns(33, 39);
  if (dialect == Dialect::ILRSB) {
    --number_of_epochs;
  }

  // Header: ## record.
  read_line();
  CHECK_EQ(columns(1, 2), "##") << location;

  // Header: +␣ records.
  read_line();
  CHECK_EQ(columns(1, 2), "+ ") << location;
  number_of_satellites = integer_columns(4, 6);

  int number_of_satellite_id_records = 0;
  while (columns(1, 2) == "+ ") {
    ++number_of_satellite_id_records;
    for (int c = 10; c <= 58; c += 3) {
      auto const full_location =
          absl::StrCat(location, " columns ", c, "-", c + 2);
      if (orbits_.size() != number_of_satellites) {
        SatelliteIdentifier id;
        if (version_ == Version::A) {
          // Satellite IDs are purely numeric (and implicitly GPS) in SP3-a.
          CHECK_EQ(column(c), ' ') << full_location;
          id.group = SatelliteGroup::GPS;
        } else {
          id.group = SatelliteGroup{column(c)};
          switch (id.group) {
            case SatelliteGroup::GPS:
            case SatelliteGroup::ГЛОНАСС:
              // GPS and ГЛОНАСС (G and R) satellite IDs are supported since
              // SP3-b.
              break;
            case SatelliteGroup::General:
            case SatelliteGroup::Galileo:
            case SatelliteGroup::北斗:
            case SatelliteGroup::準天頂衛星:
            case SatelliteGroup::IRNSS:
              CHECK_GE(version_, Version::C) << full_location;
              break;
            default:
              LOG(FATAL) << "Invalid satellite identifier " << id << ": "
                         << full_location;
          }
        }
        id.index = integer_columns(c + 1, c + 2);
        CHECK_GT(id.index, 0) << full_location;
        auto const [it, inserted] = orbits_.emplace(std::piecewise_construct,
                                                    std::forward_as_tuple(id),
                                                    std::forward_as_tuple());
        CHECK(inserted) << "Duplicate satellite identifier " << id << ": "
                        << full_location;
      } else {
        CHECK_EQ(columns(c, c + 2), "  0") << full_location;
      }
    }
    read_line();
  }
  if (number_of_satellite_id_records < 5) {
    LOG(FATAL) << u8"at least 5 +␣ records expected: " << location;
  }
  if (version_ < Version::D && number_of_satellite_id_records > 5) {
    LOG(FATAL) << u8"exactly 5 +␣ records expected in SP3-" << version_ << ": "
               << location;
  }

  // Header: ++ records.
  // Ignore the satellite accuracy exponents.
  for (int i = 0; i < number_of_satellite_id_records; ++i) {
    CHECK_EQ(columns(1, 2), "++") << location;
    read_line();
  }

  // Header: first %c record.
  std::function<Instant(std::string const&)> parse_time;
  CHECK_EQ(columns(1, 2), "%c") << location;
  if (version_ < Version::C) {
    parse_time = &ParseGPSTime;
  } else {
    auto const time_system = columns(10, 12);
    if (time_system == "GLO" || time_system == "UTC") {
      parse_time = &ParseUTC;
    } else if (time_system == "TAI") {
      parse_time = &ParseTAI;
    } else if (time_system == "BDT") {
      parse_time = &Parse北斗Time;
    } else if (time_system == "GPS" || time_system == "GAL" ||
               time_system == "TAI" || time_system == "IRN" ||
               time_system == "QZS") {
      parse_time = &ParseGPSTime;
    } else {
      LOG(FATAL) << "Unexpected time system identifier " << time_system << ": "
                 << location;
    }
  }

  // Header: second %c record.
  read_line();
  CHECK_EQ(columns(1, 2), "%c") << location;

  // Header: %f records.
  read_line();
  CHECK_EQ(columns(1, 2), "%f") << location;
  read_line();
  CHECK_EQ(columns(1, 2), "%f") << location;

  // Header: %i records.
  read_line();
  CHECK_EQ(columns(1, 2), "%i") << location;
  read_line();
  CHECK_EQ(columns(1, 2), "%i") << location;

  // Header: /* records.
  read_line();
  int number_of_comment_records = 0;
  while ((dialect == Dialect::ILRSA || dialect == Dialect::ILRSB)
             ? columns(1, 3) == "%/*"
             : columns(1, 2) == "/*") {
    ++number_of_comment_records;
    read_line();
  }
  if (number_of_comment_records < 4) {
    LOG(FATAL) << "At least 4 /* records expected: " << location;
  }
  if (version_ < Version::D && number_of_comment_records > 5) {
    LOG(FATAL) << "Exactly 4 /* records expected in SP3-"
               << version_ << ": " << location;
  }

  for (int i = 0; i < number_of_epochs; ++i) {
    // *␣ record: the epoch header record.
    CHECK_EQ(columns(1, 2), "* ") << location;
    std::string epoch_string;
    if (dialect == Dialect::ILRSB) {
      int minutes = integer_columns(17, 18);
      int hours = integer_columns(14, 15);
      if (minutes == 60) {
        minutes = 0;
        ++hours;
      }
      epoch_string = absl::StrCat(
          columns(3, 6), "-", columns(8, 9), "-", columns(11, 12),
          "T", absl::Dec(hours, absl::kZeroPad2), ":",
          absl::Dec(minutes, absl::kZeroPad2), ":", columns(20, 25));
    } else {
      // Note: the seconds field is an F11.8, spanning columns 21..31, but our
      // time parser only supports milliseconds.
      epoch_string = absl::StrCat(
          columns(4, 7), "-", columns(9, 10), "-", columns(12, 13),
          "T", columns(15, 16), ":", columns(18, 19), ":", columns(21, 26));
    }
    for (char& c : epoch_string) {
      if (c == ' ') {
        c = '0';
      }
    }
    Instant const epoch = parse_time(epoch_string);
    read_line();
    for (int i = 0; i < orbits_.size(); ++i) {
      // P record: the position and clock record.
      CHECK_EQ(column(1), 'P') << location;
      SatelliteIdentifier id;
      id.group = version_ == Version::A ? SatelliteGroup::GPS
                                        : SatelliteGroup{column(2)};
      id.index = integer_columns(3, 4);
      auto const it = orbits_.find(id);
      CHECK(it != orbits_.end()) << "Unknown satellite identifier "
                                 << id << ": " << location;
      DiscreteTrajectory<ITRS>& orbit = it->second;

      Position<ITRS> const position =
          Displacement<ITRS>({float_columns(5, 18) * Kilo(Metre),
                              float_columns(19, 32) * Kilo(Metre),
                              float_columns(33, 46) * Kilo(Metre)}) +
          ITRS::origin;
      // TODO(egg): use a difference formula to compute the velocities if they
      // are not provided.
      Velocity<ITRS> velocity({NaN<Speed>(), NaN<Speed>(), NaN<Speed>()});

      read_line();
      if (version_ >= Version::C && line.has_value() && columns(1, 2) == "EP") {
        // Ignore the optional EP record (the position and clock correlation
        // record).
        read_line();
      }

      if (has_velocities_) {
        // V record: the velocity and clock rate-of-change record.
        CHECK_EQ(column(1), 'V') << location;
        if (version_ > Version::A) {
          CHECK_EQ(SatelliteGroup{column(2)}, id.group) << location;
        }
        CHECK_EQ(integer_columns(3, 4), id.index) << location;
        velocity =
            Velocity<ITRS>({float_columns(5, 18) * (Deci(Metre) / Second),
                            float_columns(19, 32) * (Deci(Metre) / Second),
                            float_columns(33, 46) * (Deci(Metre) / Second)});
      }

      read_line();
      if (version_ >= Version::C && line.has_value() && columns(1, 2) == "EV") {
        // Ignore the optional EV record (the velocity and clock rate-of-change
        // correlation record).
        read_line();
      }

      orbit.Append(epoch, {position, velocity});
    }
  }
  if (dialect != Dialect::ILRSA) {
    CHECK_EQ(columns(1, 3), "EOF") << location;
    read_line();
  }
  CHECK(!line.has_value()) << location;
}

DiscreteTrajectory<ITRS> const& StandardProduct3::orbit(
    SatelliteIdentifier const& id) const {
  return FindOrDie(orbits_, id);
}

bool operator<(StandardProduct3::SatelliteIdentifier const& left,
               StandardProduct3::SatelliteIdentifier const& right) {
  return left.group < right.group ||
         (left.group == right.group && left.index < right.index);
}

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::Version const& version) {
  return out << std::string(static_cast<char>(version), 1);
}

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::SatelliteGroup const& group) {
  return out << std::string(static_cast<char>(group), 1);
}

std::ostream& operator<<(std::ostream& out,
                         StandardProduct3::SatelliteIdentifier const& id) {
  return out << id.group << absl::StrCat(absl::Dec(id.index, absl::kZeroPad2));
}

}  // namespace internal_standard_product_3
}  // namespace astronomy
}  // namespace principia
