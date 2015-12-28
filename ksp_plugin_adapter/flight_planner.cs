using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace principia {
namespace ksp_plugin_adapter {

class FlightPlanner : ToggleableWindowRenderer {
  public FlightPlanner(ManagerInterface manager,
                       string name): base(manager) {
    this.name = name;
    burns_ = new List<BurnEditor>();
  }

  public void RenderButton() {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    if (UnityEngine.GUILayout.Button(name + "...")) {
      show_planner_ = !show_planner_;
    }
    UnityEngine.GUI.skin = old_skin;
  }

  public string name { get; set; }

  protected override string window_title {
    get {
      return name;
    }
  }

  protected override void RenderWindow() {
    if (show_planner_) {
      window_rectangle_ = UnityEngine.GUILayout.Window(
                              id         : this.GetHashCode(),
                              screenRect : window_rectangle_,
                              func       : RenderPlanner,
                              text       : name);
    }
  }

  private void RenderPlanner(int window_id) {
    var old_skin = UnityEngine.GUI.skin;
    UnityEngine.GUI.skin = null;
    UnityEngine.GUILayout.BeginVertical();
    foreach (BurnEditor burn_editor in burns_) {
      burn_editor.Render();
    } 
    if (burns_.Count > 0) {
      if (UnityEngine.GUILayout.Button(
              "Delete",
              UnityEngine.GUILayout.ExpandWidth(true))) {
        burns_.RemoveAt(burns_.Count - 1);
      }
    }
    if (UnityEngine.GUILayout.Button("Add",
                                     UnityEngine.GUILayout.ExpandWidth(true))) {
      burns_.Add(new BurnEditor());
    }
    UnityEngine.GUILayout.EndVertical();
    UnityEngine.GUI.skin = old_skin;
  }

  private bool show_planner_;
  private UnityEngine.Rect window_rectangle_;
  private List<BurnEditor> burns_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
