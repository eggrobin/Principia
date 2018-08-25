
#pragma once

#include <algorithm>
#include <tuple>
#include <utility>

#include "quantities/tuples.hpp"

namespace principia {
namespace geometry {

// These operators live in a separate (non-internal) namespace to avoid
// polluting the entire universe in cases where they are not useful.
namespace cartesian_product {

template<typename LTuple, typename RTuple>
constexpr auto operator+(LTuple const& left, RTuple const& right);

template<typename LTuple, typename RTuple>
constexpr auto operator-(LTuple const& left, RTuple const& right);

template<typename Scalar, typename Tuple,
         typename = std::enable_if_t<quantities::is_tuple_v<Tuple>>>
constexpr auto operator*(Scalar const& left, Tuple const& right);

template<typename Tuple, typename Scalar,
         typename = std::enable_if_t<quantities::is_tuple_v<Tuple>>>
constexpr auto operator*(Tuple const& left, Scalar const& right);

template<typename Scalar, typename Tuple>
constexpr auto operator/(Tuple const& left, Scalar const& right);

}  // namespace cartesian_product
}  // namespace geometry
}  // namespace principia

#include "geometry/cartesian_product_body.hpp"
