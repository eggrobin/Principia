﻿
#pragma once

#include "geometry/hilbert.hpp"

#include "geometry/grassmann.hpp"
#include "hilbert.hpp"

import principia.quantities.elementary_functions;

namespace principia {
namespace geometry {
namespace internal_hilbert {

using namespace principia::quantities::elementary_functions;

template<typename T1, typename T2>
auto Hilbert<T1, T2,
             std::enable_if_t<
                 std::conjunction_v<is_quantity<T1>, is_quantity<T2>,
                                    std::negation<std::is_same<T1, T2>>>>>::
    InnerProduct(T1 const& t1, T2 const& t2) -> InnerProductType {
  return t1 * t2;
}

template<typename T>
auto Hilbert<T, T, std::enable_if_t<is_quantity_v<T>>>::
    InnerProduct(T const& t1, T const& t2) -> InnerProductType {
  return t1 * t2;
}

template<typename T>
auto Hilbert<T, T, std::enable_if_t<is_quantity_v<T>>>::Norm²(T const& t)
    -> Norm²Type {
  return t * t;
}

template<typename T>
auto Hilbert<T, T, std::enable_if_t<is_quantity_v<T>>>::Norm(
    T const& t) -> NormType {
  return Abs(t);
}

template<typename T1, typename T2>
auto Hilbert<T1, T2,
             std::void_t<decltype(InnerProduct(std::declval<T1>(),
                                               std::declval<T2>()))>>::
    InnerProduct(T1 const& t1, T2 const& t2) -> InnerProductType {
  // Is there a better way to avoid recursion than to put our fingers inside
  // grassmann?
  return internal_grassmann::InnerProduct(t1, t2);
}

template<typename T>
auto Hilbert<T, T,
             std::void_t<decltype(InnerProduct(std::declval<T>(),
                                               std::declval<T>()))>>::
    InnerProduct(T const& t1, T const& t2) -> InnerProductType {
  // Is there a better way to avoid recursion than to put our fingers inside
  // grassmann?
  return internal_grassmann::InnerProduct(t1, t2);
}

template<typename T>
auto Hilbert<T, T,
             std::void_t<decltype(InnerProduct(std::declval<T>(),
                                               std::declval<T>()))>>::
Norm²(T const& t) -> Norm²Type {
  return t.Norm²();
}

template<typename T>
auto Hilbert<T, T,
             std::void_t<decltype(InnerProduct(std::declval<T>(),
                                               std::declval<T>()))>>::
Norm(T const& t) -> NormType {
  return t.Norm();
}

}  // namespace internal_hilbert
}  // namespace geometry
}  // namespace principia
