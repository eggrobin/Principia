export module principia.quantities;

import <pmmintrin.h>;

import <cmath>;
import <cstdio>;
import <iostream>;
import <limits>;
import <string>;
import <type_traits>;

import "base/macros.hpp";

import "base/not_constructible.hpp";
import "base/not_null.hpp";
import "base/tags.hpp";
import "serialization/quantities.pb.h";

import principia.quantities.dimensions;
import principia.quantities.generators;

using namespace principia::quantities::dimensions;
using namespace principia::quantities::generators;

using principia::base::not_null;
using principia::base::uninitialized_t;

template<typename Left, typename Right>
using Product = typename ProductGenerator<Left, Right>::Type;
template<typename Left, typename Right>
using Quotient = typename QuotientGenerator<Left, Right>::Type;

export namespace principia::quantities {

namespace internal {
template<typename Q>
constexpr Q SIUnit();
}  // namespace internal

template<typename D>
class Quantity final {
 public:
  using Dimensions = D;

  constexpr Quantity() = default;
  explicit constexpr Quantity(uninitialized_t);

  constexpr Quantity operator+() const;
  constexpr Quantity operator-() const;
  constexpr Quantity operator+(Quantity const& right) const;
  constexpr Quantity operator-(Quantity const& right) const;

  constexpr Quantity operator*(double right) const;
  constexpr Quantity operator/(double right) const;

  Quantity& operator+=(Quantity const& right);
  Quantity& operator-=(Quantity const& right);
  Quantity& operator*=(double right);
  Quantity& operator/=(double right);

  constexpr bool operator>(Quantity const& right) const;
  constexpr bool operator<(Quantity const& right) const;
  constexpr bool operator>=(Quantity const& right) const;
  constexpr bool operator<=(Quantity const& right) const;
  constexpr bool operator==(Quantity const& right) const;
  constexpr bool operator!=(Quantity const& right) const;

  void WriteToMessage(not_null<serialization::Quantity*> message) const;
  static Quantity ReadFromMessage(serialization::Quantity const& message);

 private:
  explicit constexpr Quantity(double magnitude);
  double magnitude_ = 0;

  template<typename LDimensions, typename RDimensions>
  friend constexpr typename ProductGenerator<Quantity<LDimensions>,
                           Quantity<RDimensions>>::Type operator*(
      Quantity<LDimensions> const& left,
      Quantity<RDimensions> const& right);
  template<typename LDimensions, typename RDimensions>
  friend constexpr typename QuotientGenerator<Quantity<LDimensions>,
                            Quantity<RDimensions>>::Type operator/(
      Quantity<LDimensions> const& left,
      Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend constexpr Quantity<RDimensions> operator*(
      double left,
      Quantity<RDimensions> const& right);
  template<typename RDimensions>
  friend constexpr typename QuotientGenerator<double, Quantity<RDimensions>>::Type operator/(
      double left,
      Quantity<RDimensions> const& right);

  template<typename Q>
  friend constexpr Q internal::SIUnit();

  template<typename U>
  friend __m128d ToM128D(Quantity<U> x);
};

template<typename LDimensions, typename RDimensions>
constexpr Product<Quantity<LDimensions>, Quantity<RDimensions>>
operator*(Quantity<LDimensions> const&, Quantity<RDimensions> const&);
template<typename LDimensions, typename RDimensions>
constexpr Quotient<Quantity<LDimensions>, Quantity<RDimensions>>
operator/(Quantity<LDimensions> const&, Quantity<RDimensions> const&);
template<typename RDimensions>
constexpr Quantity<RDimensions>
operator*(double, Quantity<RDimensions> const&);
template<typename RDimensions>
constexpr Quotient<double, Quantity<RDimensions>>
operator/(double, Quantity<RDimensions> const&);

// A positive infinity of |Q|.
template<typename Q>
constexpr Q Infinity = SIUnit<Q>() * std::numeric_limits<double>::infinity();
// A quiet NaN of |Q|.
template <typename Q>
CONSTEXPR_NAN Q NaN = SIUnit<Q>() * std::numeric_limits<double>::quiet_NaN();

template<typename Q>
constexpr bool IsFinite(Q const& x);

template<typename D>
std::string Format();

std::string DebugString(
    double number,
    int precision = std::numeric_limits<double>::max_digits10);
template<typename D>
std::string DebugString(
    Quantity<D> const& quantity,
    int precision = std::numeric_limits<double>::max_digits10);

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity);

}  // namespace principia::quantities

#include "quantities/quantities_body.hpp"
