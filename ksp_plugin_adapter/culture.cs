using System.Globalization;

namespace principia {
namespace ksp_plugin_adapter {

internal static class Culture {
  static Culture() {
    // Unity/Mono is screwing with the current culture, let's get unambiguous
    // conventions from a copy of the invariant culture.
    culture = new CultureInfo(""){
        NumberFormat = {
            NumberGroupSeparator = "'", PositiveInfinitySymbol = "+∞",
            // We use U+2212 MINUS SIGN rather than U+002D HYPHEN-MINUS so that
            // the signs have the same width.  When parsing. we accept both.
            NegativeSign = "−"
        }
    };
  }

  public static readonly CultureInfo culture;
}

}  // namespace ksp_plugin_adapter
}  // namespace principia
