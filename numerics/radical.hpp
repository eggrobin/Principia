
#pragma once

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace radical {
namespace internal_radical {

using quantities::Angle;
using quantities::Product;
using quantities::Quantity;
using quantities::Quotient;
using quantities::SquareRoot;

template<typename Q>
class Radical {
 public:
  constexpr Radical(Q const& real_part) : value_{real_part} {}
  Q real_part() const;
  Q imaginary_part() const;
 private:
  constexpr Radical(bool imaginary, Q const& value)
      : imaginary_{imaginary}, value_{value} {}
  bool imaginary_;
  // If |imaginary_|, the imaginary part, otherwise, the real part.
  Q value_;

  friend constexpr Radical<double> Sqrt(double x);
};

// This is only constexpr when x = -1.  This allows |i| to be constexpr.  The
// implementation is inline to allow |i| to be defined below.
constexpr Radical<double> Sqrt(double x) {
  return x == -1 ? Radical<double>{/*imaginary=*/true, 1}
                 : std::signbit(x)
                       ? Radical<double>{/*imaginary=*/true, std::sqrt(-x)}
                       : std::sqrt(x);
}

template<typename D>
Radical<SquareRoot<Quantity<D>>> Sqrt(Quantity<D> const& x);

constexpr Radical<double> i = Sqrt(-1);

template<typename Q>
Radical<Q>& operator*=(Radical<Q>& left, double right);
template<typename Q>
Radical<Q>& operator/=(Radical<Q>& left, double right);

template<typename Q1, typename Q2>
Radical<Product<Q1, Q2>> operator*(Radical<Q1> const& left, Q2 const& right);
template<typename Q1, typename Q2>
Radical<Product<Q1, Q2>> operator*(Q1 const& left, Radical<Q2> const& right);
template<typename Q1, typename Q2>
Radical<Product<Q1, Q2>> operator*(Radical<Q1> const& left,
                                   Radical<Q2> const& right);

template<typename Q1, typename Q2>
Radical<Product<Q1, Q2>> operator/(Radical<Q1> const& left, Q2 const& right);
template<typename Q1, typename Q2>
Radical<Product<Q1, Q2>> operator/(Q1 const& left, Radical<Q2> const& right);
template<typename Q1, typename Q2>
Radical<Product<Q1, Q2>> operator/(Radical<Q1> const& left,
                                   Radical<Q2> const& right);

// Argument-dependent lookup will find the following functions when using
// |Radical| arguments.  When starting from quantities, using |radical::Sin(x)|
// etc. will convert to |Radical|.

Radical<double> Sin(Radical<Angle> const& z);
Radical<double> Cos(Radical<Angle> const& z);
Radical<double> Tan(Radical<Angle> const& z);

Radical<Angle> ArcSin(Radical<double> z);
Radical<Angle> ArcCos(Radical<double> z);
Radical<Angle> ArcTan(Radical<double> w, Radical<double> z = 1);

Radical<double> Sinh(Radical<Angle> const& z);
Radical<double> Cosh(Radical<Angle> const& z);
Radical<double> Tanh(Radical<Angle> const& z);

Radical<Angle> ArcSinh(Radical<double> z);
Radical<Angle> ArcCosh(Radical<double> z);
Radical<Angle> ArcTanh(Radical<double> z);

}  // namespace internal_radical

using internal_radical::ArcCos;
using internal_radical::ArcCosh;
using internal_radical::ArcSin;
using internal_radical::ArcSinh;
using internal_radical::ArcTan;
using internal_radical::ArcTanh;
using internal_radical::Cos;
using internal_radical::Cosh;
using internal_radical::Sin;
using internal_radical::Sinh;
using internal_radical::Sqrt;
using internal_radical::Tan;
using internal_radical::Tanh;

}  // namespace radical
}  // namespace numerics
}  // namespace principia
