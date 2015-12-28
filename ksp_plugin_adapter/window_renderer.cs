using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

internal abstract class WindowRenderer : IDisposable {
  public interface ManagerInterface {
    event Action render_windows;
  }

  public WindowRenderer(ManagerInterface manager) {
    manager_ = manager;
    manager_.render_windows += RenderWindow;
  }

  ~WindowRenderer() {
    manager_.render_windows -= RenderWindow;
  }

  public void Dispose() {
    manager_.render_windows -= RenderWindow;
    GC.SuppressFinalize(this);
  }

  protected abstract void RenderWindow();

  private ManagerInterface manager_;
}

internal abstract class ToggleableWindowRenderer : WindowRenderer {
  public ToggleableWindowRenderer(ManagerInterface manager) : base(manager) {}

  protected override void RenderWindow() {
    if (show_) {
      window_rectangle_ = UnityEngine.GUILayout.Window(
                              id         : this.GetHashCode(),
                              screenRect : window_rectangle_,
                              func       : RenderContents,
                              text       : window_title);
    }
  }

  protected abstract void RenderContents(int window_id);
  protected abstract string window_title { get; }

  protected bool show_;

  private UnityEngine.Rect window_rectangle_;
}

internal struct Controlled<T> where T : IDisposable {
  public T all {
    get {
      return all_;
    }
    set {
      if (all_ != null) {
        all_.Dispose();
      }
      all_ = value;
    }
  }

  private T all_;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
