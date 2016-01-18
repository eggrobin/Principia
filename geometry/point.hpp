﻿
#pragma once

#include <iosfwd>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/barycentre_calculator.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

namespace geometry {
template <typename Vector, typename Scalar> class BarycentreCalculator;
}  // namespace geometry
namespace serialization {
class Point;
}  // namespace serialization

using quantities::is_quantity;

namespace geometry {

// Point<Vector> is an affine space on the vector space Vector. Vector should
// be equipped with operators +, -, +=, -=, ==, !=, as well as Vector * Weight
// and Vector / Weight for any Weight used in Barycentre.
template<typename Vector>
class Point {
 public:
  Point();
  ~Point() = default;

  Vector operator-(Point const& from) const;

  Point operator+(Vector const& translation) const;
  Point operator-(Vector const& translation) const;

  Point& operator+=(Vector const& translation);
  Point& operator-=(Vector const& translation);

  bool operator==(Point const& right) const;
  bool operator!=(Point const& right) const;

  void WriteToMessage(not_null<serialization::Point*> const message) const;
  static Point ReadFromMessage(serialization::Point const& message);

 private:
  Vector coordinates_;

  template<typename V>
  friend Point<V> operator+(V const& translation, Point<V> const& point);

  template<typename V>
  friend typename std::enable_if_t<is_quantity<V>::value, bool> operator<(
      Point<V> const& left, Point<V> const& right);
  template<typename V>
  friend typename std::enable_if_t<is_quantity<V>::value, bool> operator<=(
      Point<V> const& left, Point<V> const& right);
  template<typename V>
  friend typename std::enable_if_t<is_quantity<V>::value, bool> operator>=(
      Point<V> const& left, Point<V> const& right);
  template<typename V>
  friend typename std::enable_if_t<is_quantity<V>::value, bool> operator>(
      Point<V> const& left, Point<V> const& right);

  template<typename V>
  friend std::string DebugString(Point<V> const& point);

  template<typename V, typename S>
  friend class BarycentreCalculator;
};

template<typename Vector>
Point<Vector> operator+(Vector const& translation,
                        Point<Vector> const& point);

template<typename Vector>
typename std::enable_if_t<is_quantity<Vector>::value, bool> operator<(
    Point<Vector> const& left, Point<Vector> const& right);

template<typename Vector>
typename std::enable_if_t<is_quantity<Vector>::value, bool> operator<=(
    Point<Vector> const& left, Point<Vector> const& right);

template<typename Vector>
typename std::enable_if_t<is_quantity<Vector>::value, bool> operator>=(
    Point<Vector> const& left, Point<Vector> const& right);

template<typename Vector>
typename std::enable_if_t<is_quantity<Vector>::value, bool> operator>(
    Point<Vector> const& left, Point<Vector> const& right);

template<typename Vector>
std::string DebugString(Point<Vector> const& point);

template<typename Vector>
std::ostream& operator<<(std::ostream& out, Point<Vector> const& point);

// Specialize BarycentreCalculator to make it applicable to Points.
template<typename Vector, typename Weight>
class BarycentreCalculator<Point<Vector>, Weight> {
 public:
  BarycentreCalculator() = default;
  ~BarycentreCalculator() = default;

  void Add(Point<Vector> const& point, Weight const& weight);
  Point<Vector> Get() const;

 private:
  bool empty_ = true;
  decltype(std::declval<Vector>() * std::declval<Weight>()) weighted_sum_;
  Weight weight_;
};

}  // namespace geometry
}  // namespace principia

#include "geometry/point_body.hpp"
