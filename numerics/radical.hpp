
#pragma once

#include <type_traits>

#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace radical {
namespace internal_radical {

using quantities::Angle;
using quantities::Exponentiation;
using quantities::Product;
using quantities::Quantity;
using quantities::Quotient;
using quantities::Sqrt;
using quantities::SquareRoot;

template<typename Q>
class PossiblyImaginary {
 public:
  constexpr PossiblyImaginary(Q const& real_part);  // NOLINT(runtime/explicit)
  constexpr real_part() const;
  constexpr imaginary_part() const;
 private:
  constexpr PossiblyImaginary(bool imaginary, Q const& value);
  bool imaginary_;
  // If |imaginary_|, the imaginary part, otherwise, the real part.
  Q value_;

  friend constexpr PossiblyImaginary<double> Sqrt(double x);
};

template<typename Q>
Q Abs(PossiblyImaginary<Q> const& z);

// The behaviour of |Sqrt| is correct with respect to the sign of 0, i.e.,
// Sqrt(-0) is +0 * i.
template<typename D>
PossiblyImaginary<SquareRoot<Quantity<D>>> Sqrt(Quantity<D> const& x);
// This is only constexpr when x = -1, allowing |i| to be constexpr.
constexpr PossiblyImaginary<double> Sqrt(double x);

template<int exponent,
         typename Q,
         typename = std::enable_if_t<exponent % 2 == 0>>
Exponentiation<Q, exponent> Pow(PossiblyImaginary<Q> z);

// Equality is correct with respect to the sign of 0, i.e., -0, +0, -0 * i, and
// +0 * i are equal.
template<typename Q>
bool operator==(PossiblyImaginary<Q> const& left,
                PossiblyImaginary<Q> const& right);
template<typename Q>
bool operator!=(PossiblyImaginary<Q> const& left,
                PossiblyImaginary<Q> const& right);

template<typename Q>
PossiblyImaginary<Q> operator-(PossiblyImaginary<Q> const& z);

template<typename Q1, typename Q2>
PossiblyImaginary<Product<Q1, Q2>> operator*(PossiblyImaginary<Q1> const& left,
                                             Q2 const& right);
template<typename Q1, typename Q2>
PossiblyImaginary<Product<Q1, Q2>> operator*(
    Q1 const& left,
    PossiblyImaginary<Q2> const& right);
template<typename Q1, typename Q2>
PossiblyImaginary<Product<Q1, Q2>> operator*(
    PossiblyImaginary<Q1> const& left,
    PossiblyImaginary<Q2> const& right);

template<typename Q1, typename Q2>
PossiblyImaginary<Product<Q1, Q2>> operator/(PossiblyImaginary<Q1> const& left,
                                             Q2 const& right);
template<typename Q1, typename Q2>
PossiblyImaginary<Product<Q1, Q2>> operator/(
    Q1 const& left,
    PossiblyImaginary<Q2> const& right);
template<typename Q1, typename Q2>
PossiblyImaginary<Product<Q1, Q2>> operator/(
    PossiblyImaginary<Q1> const& left,
    PossiblyImaginary<Q2> const& right);

template<typename Q>
PossiblyImaginary<Q>& operator*=(PossiblyImaginary<Q>& left, double right);
template<typename Q>
PossiblyImaginary<Q>& operator/=(PossiblyImaginary<Q>& left, double right);

// Argument-dependent lookup will find the following functions when using
// |PossiblyImaginary| arguments.  When starting from real quantities, using
// |radical::Sin(x)| etc. will convert to |PossiblyImaginary|.

PossiblyImaginary<double> Sin(PossiblyImaginary<Angle> const& z);
PossiblyImaginary<double> Cos(PossiblyImaginary<Angle> const& z);
PossiblyImaginary<double> Tan(PossiblyImaginary<Angle> const& z);

PossiblyImaginary<Angle> ArcSin(PossiblyImaginary<double> const& z);
PossiblyImaginary<Angle> ArcCos(PossiblyImaginary<double> const& z);
PossiblyImaginary<Angle> ArcTan(PossiblyImaginary<double> const& w,
                                PossiblyImaginary<double> const& z = 1);

PossiblyImaginary<double> Sinh(PossiblyImaginary<Angle> const& z);
PossiblyImaginary<double> Cosh(PossiblyImaginary<Angle> const& z);
PossiblyImaginary<double> Tanh(PossiblyImaginary<Angle> const& z);

PossiblyImaginary<Angle> ArcSinh(PossiblyImaginary<double> const& z);
PossiblyImaginary<Angle> ArcCosh(PossiblyImaginary<double> const& z);
PossiblyImaginary<Angle> ArcTanh(PossiblyImaginary<double> const& z);

// Some implementations to allow |i| to be defined constexpr below.

template<typename Q>
constexpr PossiblyImaginary<Q>::PossiblyImaginary(Q const& real_part)
    : value_{real_part} {}

template<typename Q>
constexpr PossiblyImaginary<Q>::PossiblyImaginary(bool imaginary,
                                                  Q const& value)
    : imaginary_{imaginary}, value_{value} {}

constexpr PossiblyImaginary<double> Sqrt(double x) {
  return x == -1 ? PossiblyImaginary<double>{/*imaginary=*/true, 1}
                 : std::signbit(x)
                       ? PossiblyImaginary<double>{/*imaginary=*/true, Sqrt(-x)}
                       : Sqrt(x);
}

constexpr PossiblyImaginary<double> i = Sqrt(-1);

template<typename Q>
std::string DebugString(PossiblyImaginary<Q> const& z);

template<typename Q>
std::ostream& operator<<(std::ostream& out, PossiblyImaginary<Q> const& z);

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
