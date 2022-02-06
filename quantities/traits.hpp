
#pragma once

#include <type_traits>

#include "base/not_constructible.hpp"

import principia.quantities;

namespace principia {
namespace quantities {
namespace internal_traits {

using base::not_constructible;

// A type trait for testing if a type is a quantity.
template<typename T>
struct is_quantity : std::is_arithmetic<T>, not_constructible {};
template<typename D>
struct is_quantity<Quantity<D>> : std::true_type, not_constructible {};
template<typename T>
struct is_quantity<T const> : is_quantity<T> {};

template<typename T>
constexpr bool is_quantity_v = is_quantity<T>::value;

}  // namespace internal_traits

using internal_traits::is_quantity;
using internal_traits::is_quantity_v;

}  // namespace quantities
}  // namespace principia
