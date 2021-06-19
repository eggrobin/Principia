﻿
#pragma once

#include <cstddef>
#include <cstdint>
#include <ostream>

// This file has an additional layer of namespacing for several reasons:
//   - we are not entirely sure whether we really want to expose this API;
//   - |astronomy::Time| for a class that represents the time of day is weird;
//   - |astronomy::Time| would clash with the more common |quantities::Time| in
//     some parts of |astronomy|.

namespace principia {
namespace astronomy {
namespace date_time {
namespace internal_date_time {

class Date final {
 public:
  static constexpr Date YYYYMMDD(std::int64_t digits);
  static constexpr Date YYYYwwD(std::int64_t digits);
  static constexpr Date YYYYDDD(std::int64_t digits);

  static constexpr Date Calendar(int year, int month, int day);
  static constexpr Date Ordinal(int year, int day);
  static constexpr Date Week(int year, int week, int day);

  constexpr int year() const;
  constexpr int month() const;
  constexpr int day() const;

  constexpr int ordinal() const;

  constexpr int mjd() const;

  constexpr Date next_day() const;

 private:
  constexpr Date(int year, int month, int day);

  int const year_;
  int const month_;
  int const day_;
};

class Time final {
 public:
  // Checks that this represents a valid time of day as per ISO 8601, thus
  // that the components are in the normal range, or that the object represents
  // a time in a leap second, or that it represents the end of the day.
  constexpr Time(int hour, int minute, int second, int millisecond);

  static constexpr Time hhmmss_ms(int hhmmss, int ms);

  constexpr int hour() const;
  constexpr int minute() const;
  constexpr int second() const;
  constexpr int millisecond() const;

  constexpr bool is_leap_second() const;
  // Whether |*this| is 24:00:00.
  constexpr bool is_end_of_day() const;

 private:
  int const hour_;
  int const minute_;
  int const second_;
  int const millisecond_;

  friend class TimeParser;
};

class DateTime final {
 public:
  // Checks that |time| does not represent a leap second unless |date| is the
  // last day of the month.
  constexpr DateTime(Date date, Time time);

  static constexpr DateTime BeginningOfDay(Date const& date);

  constexpr Date const& date() const;
  constexpr Time const& time() const;

  // If |time()| is 24:00:00, returns an equivalent DateTime where midnight is
  // expressed as 00:00:00 on the next day; otherwise, returns |*this|.
  constexpr DateTime normalized_end_of_day() const;

 private:

  Date const date_;
  Time const time_;

  friend constexpr DateTime operator""_DateTime(char const* str,
                                                std::size_t size);
};

class JulianDate final {
 public:
  static constexpr JulianDate JD(std::int64_t digits,
                                 std::int64_t digit_count,
                                 std::int64_t fractional_digit_count);
  static constexpr JulianDate MJD(std::int64_t digits,
                                  std::int64_t digit_count,
                                  std::int64_t fractional_digit_count);

  constexpr std::int64_t day() const;
  constexpr std::int64_t fraction_numerator() const;
  constexpr std::int64_t fraction_denominator() const;

  constexpr Date CalendarDay() const;

 private:
  constexpr JulianDate(std::int64_t day,
                       std::int64_t fraction_numerator,
                       std::int64_t fraction_denominator);

  // These numbers are relative to J2000.  |fraction_denominator| is a positive
  // power of 10.
  std::int64_t const day_;
  std::int64_t const fraction_numerator_;
  std::int64_t const fraction_denominator_;
};

constexpr bool operator==(Date const& left, Date const& right);
constexpr bool operator!=(Date const& left, Date const& right);
constexpr bool operator<(Date const& left, Date const& right);
constexpr bool operator>(Date const& left, Date const& right);
constexpr bool operator<=(Date const& left, Date const& right);
constexpr bool operator>=(Date const& left, Date const& right);
constexpr Date operator""_Date(char const* str, std::size_t size);
std::ostream& operator<<(std::ostream& out, Date const& date);

constexpr bool operator==(Time const& left, Time const& right);
constexpr bool operator!=(Time const& left, Time const& right);
constexpr Time operator""_Time(char const* str, std::size_t size);
std::ostream& operator<<(std::ostream& out, Time const& time);

constexpr bool operator==(DateTime const& left, DateTime const& right);
constexpr bool operator!=(DateTime const& left, DateTime const& right);
constexpr DateTime operator""_DateTime(char const* str, std::size_t size);
std::ostream& operator<<(std::ostream& out, DateTime const& date_time);

// Returns true if the string can be interpreted as a Julian date.
constexpr bool IsJulian(char const* str, std::size_t size);
constexpr JulianDate operator""_Julian(char const* str, std::size_t size);

}  // namespace internal_date_time

using internal_date_time::Date;
using internal_date_time::DateTime;
using internal_date_time::IsJulian;
using internal_date_time::JulianDate;
using internal_date_time::operator""_Date;
using internal_date_time::operator""_DateTime;
using internal_date_time::operator""_Julian;
using internal_date_time::operator""_Time;
using internal_date_time::Time;

}  // namespace date_time
}  // namespace astronomy
}  // namespace principia

#include "astronomy/date_time_body.hpp"
