using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

class ReferenceFrameSelector {

  Dictionary<CelestialBody, bool> expanded_;
  CelestialBody selected_celestial_;


  private void RenderSubtree(CelestialBody celestial, int depth) {
    bool leaf = celestial.orbitingBodies.Count == 0;
    bool expanded = expanded_[celestial];
    UnityEngine.GUILayout.Label("", UnityEngine.GUILayout.Width(10 * depth));
    UnityEngine.GUILayout.Button(leaf ? "" : (expanded ? "-" : "+"),
                                 UnityEngine.GUILayout.Width(10));
    if (UnityEngine.GUILayout.Toggle(selected_celestial_ == celestial,
                                     celestial.name)) {
      selected_celestial_ = celestial;
    }
    if (expanded) {
      foreach (CelestialBody child in celestial.orbitingBodies) {
        RenderSubtree(child, depth + 1);
      }
    }
  }

  public ReferenceFrameSelector() {
    expanded_ = new Dictionary<CelestialBody, bool>();
    foreach (CelestialBody celestial in FlightGlobals.Bodies) {
      expanded_.Add(celestial, false);
    }
    selected_celestial_ =
        FlightGlobals.currentMainBody ?? Planetarium.fetch.Sun;
    for (CelestialBody celestial = selected_celestial_.referenceBody;
         celestial != null;
         celestial = celestial.referenceBody) {
      expanded_[FlightGlobals.currentMainBody] = true;
    }
  }

  public void Render() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    UnityEngine.GUILayout.BeginVertical();
    RenderSubtree(celestial : Planetarium.fetch.Sun, depth : 0);
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.skin = old_skin
  }
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
