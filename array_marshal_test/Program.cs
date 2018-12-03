using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

internal struct KeplerianElements {}

class Program {

  [DllImport("kernel32", CharSet = CharSet.Unicode, SetLastError = true)]
  [return: MarshalAs(UnmanagedType.Bool)]
  static extern bool SetDllDirectory(string lpPathName);
  static void Main(string[] args) {
  
    if (!SetDllDirectory(@"C:\Users\robin\Projects\Kerbal Space Program\Plugins\Principia\Release\GameData\Principia\x64")) {
      throw new Exception( "Failed to set DLL directory (error code " +
              Marshal.GetLastWin32Error() + ").");
    }

    var plugin_ = Interface.NewPlugin(
          game_epoch:"2000-01-01T12:00:00",
          solar_system_epoch:"2000-01-01T12:00:00",
          planetarium_rotation_in_degrees: 0);
    var parameters = new Interface.BodyParameters{angular_frequency = "1 deg / s",
                                 axis_declination = "1 deg",
                                 axis_right_ascension = "1 deg",
                                 gravitational_parameter = "1 km^3/s^2",
                                 mean_radius = "1 km",
                                 name = "test_body",
                                 reference_angle = "42 deg",
                                 reference_instant = "2000-01-01T12:34:56",
                                 reference_radius="1 m"};
    parameters.geopotential = new Interface.BodyGeopotentialElement[]{
    new Interface.BodyGeopotentialElement{degree="2", order="0", cos="2",sin="0"},
    new Interface.BodyGeopotentialElement{degree="2", order="2", cos="3",sin="4"},
    new Interface.BodyGeopotentialElement{degree="3", order="0", cos="3",sin= "0" } };
        plugin_.InsertCelestialAbsoluteCartesian(
    celestial_index:0,
    parent_index:null,
    body_parameters:parameters,
    x:"1 m",
    y:"1 m",
    z:"1 m",
    vx:"1 m/s",
    vy:"1 m/s",
    vz:"1 m/s");
  }
}
internal static partial class Interface {


[StructLayout(LayoutKind.Sequential)]
internal partial struct BodyGeopotentialElement {
  public String degree;
  public String order;
  public String cos;
  public String sin;
}

[StructLayout(LayoutKind.Sequential)]
internal partial struct BodyParameters {
  public String name;
  public String gravitational_parameter;
  public String reference_instant;
  public String mean_radius;
  public String axis_right_ascension;
  public String axis_declination;
  public String reference_angle;
  public String angular_frequency;
  public String reference_radius;
  public String j2;
  [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(BodyGeopotentialElementArrayMarshaler))]
  private BodyGeopotentialElement[] geopotential_;
  private int geopotential_size_;
  public BodyGeopotentialElement[] geopotential {
    get { return geopotential_; }
    set { Console.WriteLine(value?.Length.ToString() ?? "null"); geopotential_ = value; geopotential_size_ = value?.Length ?? 0; }
  }
}

const string dll_path = "principia.dll";

  [DllImport(dllName           : dll_path,
             EntryPoint        = "principia__NewPlugin",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern IntPtr NewPlugin(
      [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(InUTF8Marshaler))] String game_epoch,
      [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(InUTF8Marshaler))] String solar_system_epoch,
      double planetarium_rotation_in_degrees);

  [DllImport(dllName           : dll_path,
             EntryPoint        = "principia__InsertCelestialAbsoluteCartesian",
             CallingConvention = CallingConvention.Cdecl)]
  internal static extern void InsertCelestialAbsoluteCartesian(
      this IntPtr plugin,
      int celestial_index,
      [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(OptionalMarshaler<int>))] BoxedInt32 parent_index,
      BodyParameters body_parameters,
      [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(InUTF8Marshaler))] String x,
      [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(InUTF8Marshaler))] String y,
      [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(InUTF8Marshaler))] String z,
      [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(InUTF8Marshaler))] String vx,
      [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(InUTF8Marshaler))] String vy,
      [MarshalAs(UnmanagedType.CustomMarshaler, MarshalTypeRef = typeof(InUTF8Marshaler))] String vz);
}

}
}
