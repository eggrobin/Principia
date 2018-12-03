using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace principia {
namespace ksp_plugin_adapter {

class Program {
  static void Main(string[] args) {
    var plugin_ = Interface.NewPlugin(
          game_epoch:"2000-01-01T12:00:00",
          solar_system_epoch:"2000-01-01T12:00:00",
          planetarium_rotation_in_degrees: 0);
    new Interface.BodyParameters{angular_frequency = "1 deg / s",
                                 axis_declination = "1 deg",
                                 axis_right_ascension = "1 deg",
                                 gravitational_parameter = "1 km^3/s^2",
                                 mean_radius = "1 km",
                                 name = "test_body",
                                 reference_angle = "42 deg",
                                 reference_instant = "2000-01-01T12:34:56"};
    plugin_.InsertCelestialAbsoluteCartesian(
    celestial_index:0,
    parent_index:null,
    body_parameters:parameters,)
  }
}

internal abstract class UTF8Marshaler : ICustomMarshaler {
  public abstract void CleanUpNativeData(IntPtr native_data);
  public abstract IntPtr MarshalManagedToNative(object managed_object);
  public abstract object MarshalNativeToManaged(IntPtr native_data);

  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  protected readonly static Encoding utf8_ =
      new UTF8Encoding(encoderShouldEmitUTF8Identifier : false,
                       throwOnInvalidBytes             : true);
}

// A marshaler for in parameter UTF-8 strings whose ownership is not taken from
// the caller.
internal class InUTF8Marshaler : UTF8Marshaler {
  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  public override void CleanUpNativeData(IntPtr native_data) {
    Marshal.FreeHGlobal(native_data);
  }

  public override IntPtr MarshalManagedToNative(object managed_object) {
    var value = managed_object as String;
    if (value == null) {
      throw new Exception(String.Format(CultureInfo.InvariantCulture,
                                    "|{0}| must be used on a |{1}|.",
                                    GetType().Name,
                                    typeof(String).Name));
    }
    int size = utf8_.GetByteCount(value);
    IntPtr buffer = Marshal.AllocHGlobal(size + 1);
    while (bytes_.Length < size + 1) {
      bytes_ = new byte[2 * bytes_.Length];
    }
    utf8_.GetBytes(value, 0, value.Length, bytes_, 0);
    bytes_[size] = 0;
    Marshal.Copy(bytes_, 0, buffer, size + 1);
    return buffer;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    throw new Exception("use |OutUTF8Marshaler| for out parameters");
  }

  private readonly static InUTF8Marshaler instance_ = new InUTF8Marshaler();
  private byte[] bytes_ = new byte[1];
}

internal static partial class Interface {

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
  public String j2;
  public String reference_radius;
}

const string dll_path = "ksp_plugin.dll";

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
