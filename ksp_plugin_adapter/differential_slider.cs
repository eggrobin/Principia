﻿using System;
using System.Globalization;

namespace principia {
namespace ksp_plugin_adapter {

internal class DifferentialSlider : ScalingRenderer {
  public delegate string ValueFormatter(double value);

  // Rates are in units of |value| per real-time second.
  public DifferentialSlider(string label,
                            string unit,
                            double log10_lower_rate,
                            double log10_upper_rate,
                            double min_value = double.NegativeInfinity,
                            double max_value = double.PositiveInfinity,
                            ValueFormatter formatter = null,
                            UnityEngine.Color? text_colour = null) {
    label_ = label;
    unit_ = unit;
    if (formatter == null) {
      format_ = v => v.ToString("#,0.000", Culture.culture);
    } else {
      format_ = formatter;
    }
    log10_lower_rate_ = log10_lower_rate;
    log10_upper_rate_ = log10_upper_rate;
    min_value_ = min_value;
    max_value_ = max_value;
    text_colour_ = text_colour;
  }

  public double value { get; set; }

  // Renders the |DifferentialSlider|.  Returns true if and only if |value|
  // changed.
  public bool Render(bool enabled) {
    bool value_changed = false;

    using (new UnityEngine.GUILayout.HorizontalScope()) {
      {
        var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.label);
        if (text_colour_.HasValue) {
          style.normal.textColor = text_colour_.Value;
        }
        UnityEngine.GUILayout.Label(text    : label_,
                                    options : GUILayoutWidth(3),
                                    style   : style);
      }

      {
        var style = new UnityEngine.GUIStyle(UnityEngine.GUI.skin.label);
        style.alignment = UnityEngine.TextAnchor.MiddleRight;
        UnityEngine.GUILayout.Label(
            text    : format_(value),
            style   : style,
            options : GUILayoutWidth(5 + (unit_ == null ? 2 : 0)));
      }
      UnityEngine.GUILayout.Label(
          text    : unit_ ?? "",
          options : GUILayoutWidth(unit_ == null ? 0 : 2));

      if (enabled) {
        if (!UnityEngine.Input.GetMouseButton(0)) {
          slider_position_ = 0;
        }
        slider_position_ = UnityEngine.GUILayout.HorizontalSlider(
            value      : slider_position_,
            leftValue  : -1,
            rightValue : 1,
            options    : UnityEngine.GUILayout.ExpandWidth(true));

        if (UnityEngine.GUILayout.Button("0", GUILayoutWidth(1))) {
          value_changed = true;
          value = 0;
        }
        if (slider_position_ != 0.0) {
          value_changed = true;
          value += Math.Sign(slider_position_) *
                   Math.Pow(10, log10_lower_rate_ +
                                    (log10_upper_rate_ - log10_lower_rate_) *
                                        Math.Abs(slider_position_)) *
                   (DateTime.Now - last_time_).TotalSeconds;
          value = Math.Min(Math.Max(min_value_, value), max_value_);
        }
      } else {
        slider_position_ = 0;
      }
      last_time_ = DateTime.Now;
    }
    return value_changed;
  }

  private float slider_position_ = 0.0f;
  private DateTime last_time_;

  private readonly string label_;
  private readonly string unit_;

  private readonly double log10_lower_rate_ = -3;
  private readonly double log10_upper_rate_ = 3.5;
  private readonly double min_value_;
  private readonly double max_value_;

  private readonly ValueFormatter format_;
  private readonly UnityEngine.Color? text_colour_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
