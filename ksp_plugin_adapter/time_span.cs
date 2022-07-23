﻿using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text.RegularExpressions;

namespace principia {
namespace ksp_plugin_adapter {

// This class starts with "Principia" to avoid confusion with the .Net TimeSpan.
class PrincipiaTimeSpan {
  public PrincipiaTimeSpan(double seconds) {
    seconds_ = seconds;
  }

  public bool Split(out int days,
                    out int hours,
                    out int minutes,
                    out double seconds) {
    days = 0;
    hours = 0;
    minutes = 0;
    seconds = seconds_;
    try {
      seconds = seconds_ % date_time_formatter.Minute;
      minutes = ((int)(seconds_ - seconds) % date_time_formatter.Hour) /
                date_time_formatter.Minute;
      hours =
          ((int)(seconds_ - seconds - minutes * date_time_formatter.Minute) %
           date_time_formatter.Day) /
          date_time_formatter.Hour;
      days = (int)(seconds_ -
                   seconds -
                   minutes * date_time_formatter.Minute -
                   hours * date_time_formatter.Hour) /
             date_time_formatter.Day;
      return true;
    } catch (OverflowException) {
      return false;
    }
  }

  // Formats a duration, optionally omitting leading components if they are 0,
  // and leading 0s on the days; optionally exclude seconds.
  public string Format(bool with_leading_zeroes, bool with_seconds) {
    return seconds_.ToString("+;−") +
           FormatPositive(with_leading_zeroes, with_seconds);
  }

  public string FormatPositive(
      bool with_leading_zeroes,
      bool with_seconds,
      int second_decimals = 1) {
    if (!Split(out int days,
               out int hours,
               out int minutes,
               out double seconds)) {
      // In case of error, saturate to the largest representable value.
      days = int.MaxValue;
      hours = 0;
      minutes = 0;
      seconds = 0;
    }
    var components = new List<string>();
    if (with_leading_zeroes) {
      components.Add(day_is_short
                         ? days.ToString("0000;0000")
                         : days.ToString("000;000"));
    } else if (days != 0) {
      components.Add(days.ToString("0;0"));
    }
    if (components.Count > 0) {
      components.Add($"{nbsp}{day_symbol}{nbsp}");
    }
    if (components.Count > 0 || with_leading_zeroes || hours != 0) {
      components.Add(day_is_short
                         ? hours.ToString("0;0")
                         : hours.ToString("00;00"));
      components.Add($"{nbsp}h{nbsp}");
    }
    if (components.Count > 0 ||
        with_leading_zeroes ||
        minutes != 0 ||
        !with_seconds) {
      components.Add(minutes.ToString("00;00"));
      components.Add($"{nbsp}min");
    }
    if (with_seconds) {
      string fractional_0s = new string('0', second_decimals);
      string seconds_format = $"00.{fractional_0s};00.{fractional_0s}";
      components.Add($"{nbsp}{seconds.ToString(seconds_format)}{nbsp }s");
    }
    return string.Join("", components.ToArray());
  }

  public double total_seconds => seconds_;

  public static bool TryParse(string text,
                              out PrincipiaTimeSpan time_span) {
    time_span = new PrincipiaTimeSpan(double.NaN);
    // Using a technology that is customarily used to parse HTML.
    string pattern = @"^[+]?\s*(?:(?<days>\d+)\s*" +
                     day_symbol +
                     @"\s*)?" +
                     @"(?:(?<hours>\d+)\s*h\s*)?" +
                     @"(?:(?<minutes>\d+)\s*min\s*)?" +
                     @"(?:(?<seconds>[0-9.,']+)\s*s\s*)?$";
    var regex = new Regex(pattern);
    var match = regex.Match(text);
    if (!match.Success) {
      return false;
    }

    var days_group = match.Groups["days"];
    var hours_group = match.Groups["hours"];
    var minutes_group = match.Groups["minutes"];
    var seconds_group = match.Groups["seconds"];
    string days = days_group.Success ? days_group.Value : "0";
    string hours = hours_group.Success ? hours_group.Value : "0";
    string minutes = minutes_group.Success ? minutes_group.Value : "0";
    string seconds = seconds_group.Success ? seconds_group.Value : "0";
    if (!int.TryParse(days, out int d) ||
        !int.TryParse(hours, out int h) ||
        !int.TryParse(minutes, out int min) ||
        !double.TryParse(seconds.Replace(',', '.'),
                         NumberStyles.AllowDecimalPoint |
                         NumberStyles.AllowThousands,
                         Culture.culture.NumberFormat,
                         out double s)) {
      return false;
    }
    time_span = new PrincipiaTimeSpan(
        Convert.ToDouble(d) * date_time_formatter.Day +
        Convert.ToDouble(h) * date_time_formatter.Hour +
        Convert.ToDouble(min) * date_time_formatter.Minute +
        s);
    return true;
  }

  public static int day_duration => date_time_formatter.Day;

  public static string day_symbol =>
      hour_divides_day && day_is_short
          ? "d" + (int)(date_time_formatter.Day / date_time_formatter.Hour)
          : "d";

  private static bool day_is_short =>
      date_time_formatter.Day / date_time_formatter.Hour < 10;

  private static bool hour_divides_day =>
      (int)(date_time_formatter.Day / date_time_formatter.Hour) *
      date_time_formatter.Hour ==
      date_time_formatter.Day;

  private static IDateTimeFormatter date_time_formatter =>
      KSPUtil.dateTimeFormatter;

  private readonly double seconds_;
  private const string nbsp = "\xA0";
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
