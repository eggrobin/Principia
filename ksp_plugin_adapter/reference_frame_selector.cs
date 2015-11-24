using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

static class CelestialExtensions {
  public static bool is_leaf(this CelestialBody celestial) {
    return celestial.orbitingBodies.Count == 0;
  }

  public static bool is_root(this CelestialBody celestial) {
    return celestial.orbit == null;
  }
}

class ReferenceFrameSelector {

  Dictionary<CelestialBody, bool> expanded_;
  CelestialBody selected_celestial_;

  private void RenderSubtree(CelestialBody celestial, int depth) {
    const int offset = 30;
    UnityEngine.GUILayout.BeginHorizontal();
    if (!celestial.is_root()) {
      UnityEngine.GUILayout.Label(
          "",
          UnityEngine.GUILayout.Width(offset * (depth - 1)));
      if (UnityEngine.GUILayout.Button(
              celestial.is_leaf() ? "" : (expanded_[celestial] ? "-" : "+"),
              UnityEngine.GUILayout.Width(offset))) {
        expanded_[celestial] = !expanded_[celestial];
      }
    }
    if (UnityEngine.GUILayout.Toggle(selected_celestial_ == celestial,
                                     celestial.name)) {
      selected_celestial_ = celestial;
    }
    UnityEngine.GUILayout.EndHorizontal();
    if (celestial.is_root() || (!celestial.is_leaf() && expanded_[celestial])) {
      foreach (CelestialBody child in celestial.orbitingBodies) {
        RenderSubtree(child, depth + 1);
      }
    }
  }

  public ReferenceFrameSelector() {
    expanded_ = new Dictionary<CelestialBody, bool>();
    foreach (CelestialBody celestial in FlightGlobals.Bodies) {
      if (!celestial.is_leaf() && !celestial.is_root()) {
        expanded_.Add(celestial, false);
      }
    }
    selected_celestial_ =
        FlightGlobals.currentMainBody ?? Planetarium.fetch.Sun;
    for (CelestialBody celestial = selected_celestial_;
         celestial.orbit != null;
         celestial = celestial.orbit.referenceBody) {
      if (!celestial.is_leaf()) {
        expanded_[celestial] = true;
      }
    }
  }

  public void Render() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    UnityEngine.GUILayout.BeginVertical();
    RenderSubtree(celestial : Planetarium.fetch.Sun, depth : 0);
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.skin = old_skin;
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
