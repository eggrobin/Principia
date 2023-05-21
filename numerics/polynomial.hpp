#pragma once

#include <algorithm>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "base/macros.hpp"
#include "base/not_null.hpp"
#include "base/traits.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/point.hpp"
#include "geometry/traits.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"
#include "serialization/numerics.pb.h"

// The presence of an operator+ below causes a bizarre compilation error in
// seemingly unrelated code in PolynomialTest.VectorSpace.
#define PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS \
  !PRINCIPIA_COMPILER_MSVC || !(_MSC_FULL_VER == 192'930'036 || \
                                _MSC_FULL_VER == 192'930'037 || \
                                _MSC_FULL_VER == 192'930'038 || \
                                _MSC_FULL_VER == 192'930'133 || \
                                _MSC_FULL_VER == 192'930'139 || \
                                _MSC_FULL_VER == 192'930'143 || \
                                _MSC_FULL_VER == 192'930'147 || \
                                _MSC_FULL_VER == 193'431'937 || \
                                _MSC_FULL_VER == 193'431'942 || \
                                _MSC_FULL_VER == 193'431'944 || \
                                _MSC_FULL_VER == 193'532'216 || \
                                _MSC_FULL_VER == 193'532'217 || \
                                _MSC_FULL_VER == 193'632'532)

namespace principia {
namespace numerics {
FORWARD_DECLARE_FROM(
    polynomial,
    TEMPLATE(typename Value, typename Argument, int degree_,
            template<typename, typename, int> typename Evaluator) class,
    PolynomialInMonomialBasis);
}  // namespace numerics

namespace mathematica {
FORWARD_DECLARE_FUNCTION_FROM(
    mathematica,
    TEMPLATE(typename Value, typename Argument, int degree_,
             template<typename, typename, int> typename Evaluator,
             typename OptionalExpressIn) std::string,
    ToMathematicaBody,
    (numerics::_polynomial::
         PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
             polynomial,
     OptionalExpressIn express_in));
}  // namespace mathematica

namespace numerics {
namespace _polynomial {
namespace internal {

using namespace principia::base::_not_constructible;
using namespace principia::base::_not_null;
using namespace principia::base::_traits;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_point;
using namespace principia::geometry::_traits;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_tuples;

// |Value_| must belong to an affine space.  |Argument_| must belong to a ring
// or to Point based on a ring.
template<typename Value_, typename Argument_>
class Polynomial {
 public:
  using Argument = Argument_;
  using Value = Value_;

  // This virtual destructor makes this class and its subclasses non-literal, so
  // constexpr-ness is a bit of a lie for polynomials.
  // TODO(phl): Consider providing an explicit deleter function that would allow
  // making the destructor protected and nonvirtual.
  virtual ~Polynomial() = default;

  virtual Value operator()(Argument const& argument) const = 0;
  virtual Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const = 0;

  // Only useful for benchmarking, analyzing performance or for downcasting.  Do
  // not use in other circumstances.
  virtual int degree() const = 0;

  // Only useful for logging.  Do not use in real code.
  virtual bool is_zero() const = 0;

  virtual void WriteToMessage(
      not_null<serialization::Polynomial*> message) const = 0;

  // The evaluator is not part of the serialization because it's fine to read
  // with a different evaluator than the one the polynomial was written with.
  template<template<typename, typename, int> typename Evaluator>
  static not_null<std::unique_ptr<Polynomial>> ReadFromMessage(
      serialization::Polynomial const& message);
};

template<typename Value_, typename Argument_, int degree_,
         template<typename, typename, int> typename Evaluator>
class PolynomialInMonomialBasis : public Polynomial<Value_, Argument_> {
 public:
  using Argument = Argument_;
  using Value = Value_;

  // Equivalent to:
  //   std::tuple<Value,
  //              Derivative<Value, Argument>,
  //              Derivative<Derivative<Value, Argument>>...>
  using Coefficients = Derivatives<Value, Argument, degree_ + 1>;

  // The coefficients are relative to origin; in other words they are applied to
  // powers of (argument - origin).
  constexpr PolynomialInMonomialBasis(Coefficients coefficients,
                                      Argument const& origin);
  template<typename A = Argument,
           typename = std::enable_if_t<is_vector_v<A>>>
  explicit constexpr PolynomialInMonomialBasis(
      Coefficients coefficients);

  // A polynomial may be explicitly converted to a higher degree (possibly with
  // a different evaluator).
  template<int higher_degree_,
           template<typename, typename, int> class HigherEvaluator>
  explicit operator PolynomialInMonomialBasis<
      Value, Argument, higher_degree_, HigherEvaluator>() const;

  Value operator()(Argument const& argument) const override;
  Derivative<Value, Argument> EvaluateDerivative(
      Argument const& argument) const override;

  constexpr int degree() const override;
  bool is_zero() const override;

  Argument const& origin() const;

  // Returns a copy of this polynomial adjusted to the given origin.
  PolynomialInMonomialBasis AtOrigin(Argument const& origin) const;

  template<int order = 1>
  PolynomialInMonomialBasis<
      Derivative<Value, Argument, order>, Argument, degree_ - order, Evaluator>
  Derivative() const;

  // The constant term of the result is zero.
  template<typename V = Value,
           typename = std::enable_if_t<is_vector_v<V>>>
  PolynomialInMonomialBasis<Primitive<Value, Argument>,
                            Argument, degree_ + 1, Evaluator>
  Primitive() const;

  template<typename V = Value,
           typename = std::enable_if_t<is_vector_v<V>>>
  quantities::_named_quantities::Primitive<Value, Argument> Integrate(
      Argument const& argument1,
      Argument const& argument2) const;

  PolynomialInMonomialBasis& operator+=(const PolynomialInMonomialBasis& right);
  PolynomialInMonomialBasis& operator-=(const PolynomialInMonomialBasis& right);

  void WriteToMessage(
      not_null<serialization::Polynomial*> message) const override;
  static PolynomialInMonomialBasis ReadFromMessage(
      serialization::Polynomial const& message);

 private:
  Coefficients coefficients_;
  Argument origin_;

  template<typename V, typename A, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, r, E>
  friend operator-(PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename V, typename A, int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, PRINCIPIA_MAX(l, r), E>
  friend operator+(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename V, typename A, int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, PRINCIPIA_MAX(l, r), E>
  friend operator-(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename S,
           typename V, typename A, int d,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Product<S, V>, A, d, E>
  friend operator*(S const& left,
                   PolynomialInMonomialBasis<V, A, d, E> const& right);
  template<typename S,
           typename V, typename A, int d,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Product<V, S>, A, d, E>
  friend operator*(PolynomialInMonomialBasis<V, A, d, E> const& left,
                   S const& right);
  template<typename S,
           typename V, typename A, int d,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Quotient<V, S>, A, d, E>
  friend operator/(PolynomialInMonomialBasis<V, A, d, E> const& left,
                   S const& right);
  template<typename L, typename R, typename A,
           int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Product<L, R>, A, l + r, E>
  friend operator*(
      PolynomialInMonomialBasis<L, A, l, E> const& left,
      PolynomialInMonomialBasis<R, A, r, E> const& right);
#if PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS
  template<typename V, typename A, int l,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, l, E>
  friend operator+(
      PolynomialInMonomialBasis<Difference<V>, A, l, E> const& left,
      V const& right);
#endif
  template<typename V, typename A, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<V, A, r, E>
  friend operator+(
      V const& left,
      PolynomialInMonomialBasis<Difference<V>, A, r, E> const& right);
  template<typename V, typename A, int l,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Difference<V>, A, l, E>
  friend operator-(PolynomialInMonomialBasis<V, A, l, E> const& left,
                   V const& right);
  template<typename V, typename A, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<Difference<V>, A, r, E>
  friend operator-(V const& left,
                   PolynomialInMonomialBasis<V, A, r, E> const& right);
  template<typename L, typename R, typename A,
           int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<L, A, l * r, E>
  friend Compose(PolynomialInMonomialBasis<L, R, l, E> const& left,
                 PolynomialInMonomialBasis<R, A, r, E> const& right);
  template<typename L, typename R, typename A,
           int l, int r,
           template<typename, typename, int> typename E>
  constexpr PolynomialInMonomialBasis<
      typename Hilbert<L, R>::InnerProductType, A, l + r, E>
  friend PointwiseInnerProduct(
      PolynomialInMonomialBasis<L, A, l, E> const& left,
      PolynomialInMonomialBasis<R, A, r, E> const& right);
  template<typename V, typename A, int d,
           template<typename, typename, int> typename E>
  friend std::ostream& operator<<(
      std::ostream& out,
      PolynomialInMonomialBasis<V, A, d, E> const& polynomial);
  template<typename V, typename A, int d,
           template<typename, typename, int> class E,
           typename O>
  friend std::string mathematica::_mathematica::internal::ToMathematicaBody(
      PolynomialInMonomialBasis<V, A, d, E> const& polynomial,
      O express_in);
};

// Vector space of polynomials.

template<typename Value, typename Argument, int rdegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator>
operator+(PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
              right);

template<typename Value, typename Argument, int rdegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator>
operator-(PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
              right);

template<typename Value, typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument,
                                    PRINCIPIA_MAX(ldegree_, rdegree_),
                                    Evaluator>
operator+(
    PolynomialInMonomialBasis<Value, Argument, ldegree_, Evaluator> const& left,
    PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
        right);

template<typename Value, typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument,
                                    PRINCIPIA_MAX(ldegree_, rdegree_),
                                    Evaluator>
operator-(
    PolynomialInMonomialBasis<Value, Argument, ldegree_, Evaluator> const& left,
    PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator> const&
        right);

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Product<Scalar, Value>, Argument,
                                    degree_, Evaluator>
operator*(Scalar const& left,
          PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              right);

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Product<Value, Scalar>, Argument,
                                    degree_, Evaluator>
operator*(PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              left,
          Scalar const& right);

template<typename Scalar,
         typename Value, typename Argument, int degree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Quotient<Value, Scalar>, Argument,
                                    degree_, Evaluator>
operator/(PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
              left,
          Scalar const& right);

// Algebra of polynomials.

template<typename LValue, typename RValue,
         typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Product<LValue, RValue>, Argument,
                                    ldegree_ + rdegree_, Evaluator>
operator*(
    PolynomialInMonomialBasis<LValue, Argument, ldegree_, Evaluator> const&
        left,
    PolynomialInMonomialBasis<RValue, Argument, rdegree_, Evaluator> const&
        right);

// Additive operators polynomial ± constant.

#if PRINCIPIA_COMPILER_MSVC_HANDLES_POLYNOMIAL_OPERATORS
template<typename Value, typename Argument, int ldegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, ldegree_, Evaluator>
operator+(PolynomialInMonomialBasis<Difference<Value>, Argument,
                                    ldegree_, Evaluator> const& left,
          Value const& right);
#endif

template<typename Value, typename Argument, int rdegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Value, Argument, rdegree_, Evaluator>
operator+(Value const& left,
          PolynomialInMonomialBasis<Difference<Value>, Argument,
                                    rdegree_, Evaluator> const& right);

template<typename Value, typename Argument, int ldegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Difference<Value>, Argument,
                                    ldegree_, Evaluator>
operator-(PolynomialInMonomialBasis<Value, Argument,
                                    ldegree_, Evaluator> const& left,
          Value const& right);

template<typename Value, typename Argument, int rdegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<Difference<Value>, Argument,
                                    rdegree_, Evaluator>
operator-(Value const& left,
          PolynomialInMonomialBasis<Value, Argument,
                                    rdegree_, Evaluator> const& right);

// Application monoid.

template<typename LValue, typename RValue,
         typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<LValue, Argument,
                                    ldegree_ * rdegree_, Evaluator>
Compose(
    PolynomialInMonomialBasis<LValue, RValue, ldegree_, Evaluator> const&
        left,
    PolynomialInMonomialBasis<RValue, Argument, rdegree_, Evaluator> const&
        right);

// Returns a scalar polynomial obtained by pointwise inner product of two
// vector-valued polynomials.
template<typename LValue, typename RValue,
         typename Argument, int ldegree_, int rdegree_,
         template<typename, typename, int> typename Evaluator>
constexpr PolynomialInMonomialBasis<
    typename Hilbert<LValue, RValue>::InnerProductType, Argument,
    ldegree_ + rdegree_, Evaluator>
PointwiseInnerProduct(
    PolynomialInMonomialBasis<LValue, Argument, ldegree_, Evaluator> const&
        left,
    PolynomialInMonomialBasis<RValue, Argument, rdegree_, Evaluator> const&
        right);

// Output.

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> typename Evaluator>
std::ostream& operator<<(
    std::ostream& out,
    PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator> const&
        polynomial);

}  // namespace internal

using internal::Polynomial;
using internal::PolynomialInMonomialBasis;

}  // namespace _polynomial
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_body.hpp"
