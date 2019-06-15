﻿using System;
using System.Runtime.InteropServices;

namespace principia {
namespace ksp_plugin_adapter {

internal class OutOwnedUTF16Marshaler : ICustomMarshaler {
  void ICustomMarshaler.CleanUpManagedData(object managed_object) {}

  int ICustomMarshaler.GetNativeDataSize() {
    return -1;
  }

  public static ICustomMarshaler GetInstance(string s) {
    return instance_;
  }

  void ICustomMarshaler.CleanUpNativeData(IntPtr native_data) {
    Interface.DeleteU16String(ref native_data);
  }
  
  IntPtr ICustomMarshaler.MarshalManagedToNative(object managed_object) {
    throw Log.Fatal("use |MarshalAs(UnmanagedType.LPWStr)| for in parameters");
  }

  object ICustomMarshaler.MarshalNativeToManaged(IntPtr native_data) {
    return Marshal.PtrToStringUni(native_data);
  }

  private static readonly OutOwnedUTF16Marshaler instance_ =
      new OutOwnedUTF16Marshaler();
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
