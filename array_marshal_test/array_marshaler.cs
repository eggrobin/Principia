using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {
internal static partial class Interface {
internal class BodyParametersMarshaler : ICustomMarshaler {

  [StructLayout(LayoutKind.Sequential)]
  internal partial class BodyParametersRepresentation {
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
    public IntPtr geopotential;
    public int geopotential_size;
  }

  // In addition to implementing the |ICustomMarshaler| interface, custom
  // marshalers must implement a static method called |GetInstance| that accepts
  // a |String| as a parameter and has a return type of |ICustomMarshaler|,
  // see https://goo.gl/wwmBTa.
  public static ICustomMarshaler GetInstance(String s) {
    return instance_;
  }

  public void CleanUpNativeData(IntPtr native_data) {
    Console.WriteLine("cleanup");
    var representation = new BodyParametersRepresentation();
    Marshal.PtrToStructure(native_data, representation);
    for (int i = 0; i < representation.geopotential_size; ++i) {
      Marshal.DestroyStructure<BodyGeopotentialElement>(IntPtr.Add(representation.geopotential, i * Marshal.SizeOf<BodyGeopotentialElement>()));
    }
    Marshal.FreeHGlobal(representation.geopotential);
    Marshal.FreeHGlobal(native_data);
    Console.WriteLine("Deallocated");
  }

  public IntPtr MarshalManagedToNative(object managed_object) {
    var parameters = managed_object as BodyParameters;
    var representation = new BodyParametersRepresentation {
        angular_frequency=parameters.angular_frequency,
        axis_declination=parameters.axis_declination,
        axis_right_ascension=parameters.axis_right_ascension,
        gravitational_parameter=parameters.gravitational_parameter,
        j2=parameters.j2,
        mean_radius=parameters.mean_radius,
        name=parameters.name,
        reference_angle=parameters.reference_angle,
        reference_instant=parameters.reference_instant,
        reference_radius=parameters.reference_radius};
    representation.geopotential_size = parameters.geopotential.Length;
    if (parameters.geopotential.Length == 0) {
      representation.geopotential = IntPtr.Zero;
    } else {
      int sizeof_element = Marshal.SizeOf<BodyGeopotentialElement>();
      representation.geopotential = Marshal.AllocHGlobal(sizeof_element * parameters.geopotential.Length);
      for (int i = 0; i < parameters.geopotential.Length; ++i) {
        Marshal.StructureToPtr(parameters.geopotential[i], IntPtr.Add(representation.geopotential, i * sizeof_element), false);
      }
    }
    IntPtr buffer = Marshal.AllocHGlobal(Marshal.SizeOf<BodyParametersRepresentation>());
    Marshal.StructureToPtr(representation, buffer, false);
    Console.WriteLine("Allocated " + buffer.ToInt64().ToString("X"));
    return buffer;
  }

  public object MarshalNativeToManaged(IntPtr native_data) {
    throw new Exception("use |OutUTF8Marshaler| for out parameters");
  }

  public void CleanUpManagedData(object managed_data) {
    throw new Exception("use |OutUTF8Marshaler| for out parameters");
  }
  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  private readonly static BodyParametersMarshaler instance_ = new BodyParametersMarshaler();
}

}

}  // namespace ksp_plugin_adapter
}  // namespace principia
