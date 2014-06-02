#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

#include "quantities/dimensionless.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"


namespace principia {
namespace benchmarks {

static std::size_t const dimension = 100;

inline void DimensionfulDiscreteCosineTransform(
    std::vector<quantities::Momentum>* result) {
  using quantities::Cos;
  using quantities::Dimensionless;
  using quantities::Momentum;
  using si::Radian;
  std::vector<Momentum> input(dimension);
  for (std::size_t i = 0; i < dimension; ++i) {
    input[i] = i * Momentum::SIUnit();
  }
  result->resize(dimension);
  Dimensionless sign = 1;
  Momentum sum;
  for (std::size_t k = 0; k < dimension; ++k, sign *= -1) {
    sum = Momentum();
    for (std::size_t n = 1; n < dimension - 1; ++n) {
      sum += input[n] * quantities::Cos(π * Radian / (dimension - 1) * n * k);
    }
    // We omit adding sum back to see whether the compiler will properly
    // optimise away the above loop.
    (*result)[k] = 0.5 * (input[0] + sign * input[dimension - 1]);  // + sum;
  }
}

inline void DoubleDiscreteCosineTransform(
    std::vector<double>* result) {;
  std::vector<double> input(dimension);
  for (std::size_t i = 0; i < dimension; ++i) {
    input[i] = i;
  }
  result->resize(dimension);
  double sign = 1;
  double sum;
  for (std::size_t k = 0; k < dimension; ++k, sign *= -1) {
    sum = 0;
    for (std::size_t n = 1; n < dimension - 1; ++n) {
      sum += input[n] * std::cos(M_PI / (dimension - 1) * n * k);
    }
    // We omit adding sum back to see whether the compiler will properly
    // optimise away the above loop.
    (*result)[k] = 0.5 * (input[0] + sign * input[dimension - 1]);  // + sum;
  }
}

class TestType {
 public:
  __declspec(noalias)
  TestType() = default;
  __declspec(noalias)
  TestType(double value) : value_(value) {}
  __declspec(noalias)
  double value() const { return value_; }
 private:
  double value_;
};

TestType const Pi {M_PI};

__declspec(noalias)
inline TestType Cos(TestType const& x) {
  return TestType {std::cos(x.value())};
}

__declspec(noalias)
inline TestType operator*(TestType const& x, TestType const& y) {
  return TestType {x.value() * y.value()};
}

__declspec(noalias)
inline TestType operator/(TestType const& x, TestType const& y) {
  return TestType {x.value() / y.value()};
}

__declspec(noalias)
inline TestType operator+(TestType const& x, TestType const& y) {
  return TestType {x.value() + y.value()};
}

__declspec(noalias)
inline TestType& operator*=(TestType& x, TestType const& y) {
  return x = x * y;
}

__declspec(noalias)
inline TestType& operator+=(TestType& x, TestType const& y) {
  return x = x + y;
}

inline void TestTypeDiscreteCosineTransform(
    std::vector<TestType>* result) {;
  std::vector<TestType> input(dimension);
  for (std::size_t i = 0; i < dimension; ++i) {
    input[i] = TestType {(double)i};
  }
  result->resize(dimension);
  TestType sign {1};
  TestType sum;
  for (std::size_t k = 0; k < dimension; ++k, sign *= TestType{-1}) {
    sum = TestType {0};
    for (std::size_t n = 1; n < dimension - 1; ++n) {
      sum += input[n] * Cos(Pi / TestType{(double)(dimension - 1) * n * k});
    }
    // We omit adding sum back to see whether the compiler will properly
    // optimise away the above loop.
    (*result)[k] = TestType{0.5} * (input[0] + sign * input[dimension - 1]);  // + sum;
  }
}


}  // namespace benchmarks
}  // namespace principia
