﻿
#pragma once

#include "numerics/fixed_arrays.hpp"

#include <algorithm>
#include <vector>

#include "glog/logging.h"

namespace principia {
namespace numerics {
namespace internal_fixed_arrays {

// A helper class to compute the dot product of two arrays.  |ScalarLeft| and
// |ScalarRight| are the types of the elements of the arrays.  |Left| and
// |Right| are the (deduced) types of the arrays.  They must both have an
// operator[].  |size| is the size of the arrays.
template<typename ScalarLeft, typename ScalarRight, int size, int i = size - 1>
struct DotProduct {
  template<typename Left, typename Right>
  static Product<ScalarLeft, ScalarRight> Compute(Left const& left,
                                                  Right const& right);
};

template<typename ScalarLeft, typename ScalarRight, int size>
struct DotProduct<ScalarLeft, ScalarRight, size, 0> {
  template<typename Left, typename Right>
  static Product<ScalarLeft, ScalarRight> Compute(Left const& left,
                                                  Right const& right);
};

template<typename ScalarLeft, typename ScalarRight, int size, int i>
template<typename Left, typename Right>
Product<ScalarLeft, ScalarRight>
DotProduct<ScalarLeft, ScalarRight, size, i>::Compute(Left const& left,
                                                      Right const& right) {
  return left[i] * right[i] +
         DotProduct<ScalarLeft, ScalarRight, size, i - 1>::Compute(left, right);
}

template<typename ScalarLeft, typename ScalarRight, int size>
template<typename Left, typename Right>
Product<ScalarLeft, ScalarRight>
DotProduct<ScalarLeft, ScalarRight, size, 0>::Compute(Left const& left,
                                                      Right const& right) {
  return left[0] * right[0];
}

// The |data_| member is aggregate-initialized with an empty list initializer,
// which performs value initialization on the components.  For quantities this
// calls the default constructor, for non-class types this does
// zero-initialization.
template<typename Scalar, int size_>
constexpr FixedVector<Scalar, size_>::FixedVector() : data_{} {}

template<typename Scalar, int size_>
FixedVector<Scalar, size_>::FixedVector(uninitialized_t) {}

template<typename Scalar, int size_>
constexpr FixedVector<Scalar, size_>::FixedVector(
    std::array<Scalar, size_> const& data)
    : data_(data) {}

template<typename Scalar, int size_>
constexpr FixedVector<Scalar, size_>::FixedVector(
    std::array<Scalar, size_>&& data)
    : data_(std::move(data)) {}

template<typename Scalar, int size_>
bool FixedVector<Scalar, size_>::operator==(FixedVector const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int size_>
constexpr Scalar& FixedVector<Scalar, size_>::operator[](int const index) {
  return data_[index];
}

template<typename Scalar, int size_>
constexpr Scalar const& FixedVector<Scalar, size_>::operator[](
    int const index) const {
  return data_[index];
}

template<typename Scalar, int size_>
FixedVector<Scalar, size_>::operator std::vector<Scalar>() const {
  std::vector<Scalar> result(data_.size());
  std::copy(data_.begin(), data_.end(), result.begin());
  return result;
}

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns>::FixedMatrix()
    : data_{} {}

template<typename Scalar, int rows, int columns>
FixedMatrix<Scalar, rows, columns>::FixedMatrix(uninitialized_t) {}

template<typename Scalar, int rows, int columns>
constexpr FixedMatrix<Scalar, rows, columns>::FixedMatrix(
    std::array<Scalar, rows * columns> const& data)
    : data_(data) {}

template<typename Scalar, int rows, int columns>
bool FixedMatrix<Scalar, rows, columns>::operator==(
    FixedMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows, int columns>
Scalar* FixedMatrix<Scalar, rows, columns>::operator[](int const index) {
  return &data_[index * columns];
}

template<typename Scalar, int rows, int columns>
constexpr Scalar const* FixedMatrix<Scalar, rows, columns>::operator[](
    int const index) const {
  return &data_[index * columns];
}

template<typename Scalar, int rows, int columns>
template<int r>
FixedMatrix<Scalar, rows, columns>::Row<r>::Row(const FixedMatrix* const matrix)
    : matrix_(matrix) {}

template<typename Scalar, int rows, int columns>
template<int r>
constexpr Scalar const& FixedMatrix<Scalar, rows, columns>::Row<r>::operator[](
    int index) const {
  return (matrix_->data_)[r * columns + index];
}

template<typename Scalar, int rows, int columns>
template<int r>
template<typename S>
Product<Scalar, S>
FixedMatrix<Scalar, rows, columns>::Row<r>::operator*(
    FixedVector<S, columns> const& right) {
  return DotProduct<Scalar, S, columns>::Compute(*this, right);
}

template<typename Scalar, int rows, int columns>
template<int r>
typename FixedMatrix<Scalar, rows, columns>::template Row<r>
FixedMatrix<Scalar, rows, columns>::row() const {
  return Row<r>(this);
}

template<typename ScalarLeft, typename ScalarRight, int size>
constexpr FixedVector<Difference<ScalarLeft, ScalarRight>, size> operator-(
    FixedVector<ScalarLeft, size> const& left,
    FixedVector<ScalarRight, size> const& right) {
  std::array<Difference<ScalarLeft, ScalarRight>, size> result{};
  for (int i = 0; i < size; ++i) {
    result[i] = left[i] - right[i];
  }
  return FixedVector<Difference<ScalarLeft, ScalarRight>, size>(
      std::move(result));
}

template<typename ScalarLeft, typename ScalarRight, int rows, int columns>
FixedVector<Product<ScalarLeft, ScalarRight>, rows> operator*(
    FixedMatrix<ScalarLeft, rows, columns> const& left,
    FixedVector<ScalarRight, columns> const& right) {
  std::array<Product<ScalarLeft, ScalarRight>, rows> result;
  auto const* row = left.data_.data();
  for (int i = 0; i < rows; ++i) {
    result[i] =
        DotProduct<ScalarLeft, ScalarRight, columns>::Compute(row, right.data_);
    row += columns;
  }
  return FixedVector<Product<ScalarLeft, ScalarRight>, rows>(std::move(result));
}

template<typename Scalar, int rows>
constexpr FixedStrictlyLowerTriangularMatrix<Scalar, rows>::
    FixedStrictlyLowerTriangularMatrix()
    : data_{} {}

template<typename Scalar, int rows>
FixedStrictlyLowerTriangularMatrix<Scalar, rows>::
    FixedStrictlyLowerTriangularMatrix(uninitialized_t) {}

template<typename Scalar, int rows>
constexpr FixedStrictlyLowerTriangularMatrix<Scalar, rows>::
    FixedStrictlyLowerTriangularMatrix(
        std::array<Scalar, dimension> const& data)
    : data_(data) {}

template<typename Scalar, int rows>
bool FixedStrictlyLowerTriangularMatrix<Scalar, rows>::operator==(
    FixedStrictlyLowerTriangularMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows>
Scalar* FixedStrictlyLowerTriangularMatrix<Scalar, rows>::operator[](
    int const index) {
  return &data_[index * (index - 1) / 2];
}

template<typename Scalar, int rows>
constexpr Scalar const*
FixedStrictlyLowerTriangularMatrix<Scalar, rows>::operator[](
    int const index) const {
  return &data_[index * (index - 1) / 2];
}

template<typename Scalar, int rows>
constexpr int FixedStrictlyLowerTriangularMatrix<Scalar, rows>::dimension;

template<typename Scalar, int rows>
constexpr FixedLowerTriangularMatrix<Scalar, rows>::FixedLowerTriangularMatrix()
    : data_{} {}

template<typename Scalar, int rows>
FixedLowerTriangularMatrix<Scalar, rows>::FixedLowerTriangularMatrix(
    uninitialized_t) {}

template<typename Scalar, int rows>
constexpr FixedLowerTriangularMatrix<Scalar, rows>::
    FixedLowerTriangularMatrix(std::array<Scalar, dimension> const& data)
    : data_(data) {}

template<typename Scalar, int rows>
bool FixedLowerTriangularMatrix<Scalar, rows>::operator==(
    FixedLowerTriangularMatrix const& right) const {
  return data_ == right.data_;
}

template<typename Scalar, int rows>
Scalar* FixedLowerTriangularMatrix<Scalar, rows>::operator[](
    int const index) {
  return &data_[index * (index + 1) / 2];
}

template<typename Scalar, int rows>
constexpr Scalar const*
FixedLowerTriangularMatrix<Scalar, rows>::operator[](
    int const index) const {
  return &data_[index * (index + 1) / 2];
}

template<typename Scalar, int rows>
constexpr int FixedLowerTriangularMatrix<Scalar, rows>::dimension;

}  // namespace internal_fixed_arrays
}  // namespace numerics
}  // namespace principia
