#pragma once

namespace principia {
namespace clr_benchmarks_adapter {

public ref class QuantitiesCLRBenchmark abstract sealed {
 public:
  static void DimensionfulDiscreteCosineTransform();
  static void DoubleDiscreteCosineTransform();
  static void TestTypeDiscreteCosineTransform();
};

}  // namespace clr_benchmarks_adapter
}  // namespace principia
