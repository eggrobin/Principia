using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace principia {
namespace ksp_plugin_adapter {

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
    Console.WriteLine("   Deallocated string");
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
    Console.WriteLine("   Allocated string");
    return buffer;
  }

  public override object MarshalNativeToManaged(IntPtr native_data) {
    throw new Exception("use |OutUTF8Marshaler| for out parameters");
  }

  private readonly static InUTF8Marshaler instance_ = new InUTF8Marshaler();
  private byte[] bytes_ = new byte[1];
}


}  // namespace ksp_plugin_adapter
}  // namespace principia
