﻿
#include <ostream>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/serialization.hpp"
#include "glog/logging.h"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/fixed_arrays.hpp"
#include "numerics/newhall.mathematica.h"
#include "numerics/чебышёв_series.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/numerics.pb.h"

namespace principia {

namespace geometry {
template <typename Scalar, typename Frame, int rank> class Multivector;
}  // namespace geometry
namespace numerics {
namespace internal {
template <typename Vector> class EvaluationHelper;
}  // namespace internal
}  // namespace numerics

using geometry::DoubleOrQuantityOrMultivectorSerializer;
using geometry::Multivector;
using geometry::R3Element;

namespace numerics {
namespace internal {

// The compiler does a much better job on an |R3Element<double>| than on a
// |Vector<Quantity>| so we specialize this case.
template<typename Scalar, typename Frame, int rank>
class EvaluationHelper<Multivector<Scalar, Frame, rank>> {
 public:
  EvaluationHelper(
      std::vector<Multivector<Scalar, Frame, rank>> const& coefficients,
      int const degree);
  EvaluationHelper(EvaluationHelper&& other) = default;
  EvaluationHelper& operator=(EvaluationHelper&& other) = default;

  Multivector<Scalar, Frame, rank> EvaluateImplementation(
      double const scaled_t) const;

  Multivector<Scalar, Frame, rank> coefficients(int const index) const;
  int degree() const;

 private:
  std::vector<R3Element<double>> coefficients_;
  int degree_;
};

template<typename Vector>
EvaluationHelper<Vector>::EvaluationHelper(
    std::vector<Vector> const& coefficients,
    int const degree) : coefficients_(coefficients), degree_(degree) {}

template<typename Vector>
Vector EvaluationHelper<Vector>::EvaluateImplementation(
    double const scaled_t) const {
  double const two_scaled_t = scaled_t + scaled_t;
  Vector const c_0 = coefficients_[0];
  switch (degree_) {
    case 0:
      return c_0;
    case 1:
      return c_0 + scaled_t * coefficients_[1];
    default:
      // b_degree   = c_degree.
      Vector b_i = coefficients_[degree_];
      // b_degree-1 = c_degree-1 + 2 t b_degree.
      Vector b_j = coefficients_[degree_ - 1] + two_scaled_t * b_i;
      int k = degree_ - 3;
      for (; k >= 1; k -= 2) {
        // b_k+1 = c_k+1 + 2 t b_k+2 - b_k+3.
        b_i = coefficients_[k + 1] + two_scaled_t * b_j - b_i;
        // b_k   = c_k   + 2 t b_k+1 - b_k+2.
        b_j = coefficients_[k] + two_scaled_t * b_i - b_j;
      }
      if (k == 0) {
        // b_1 = c_1 + 2 t b_2 - b_3.
        b_i = coefficients_[1] + two_scaled_t * b_j - b_i;
        // c_0 + t b_1 - b_2.
        return c_0 + scaled_t * b_i - b_j;
      } else {
        // c_0 + t b_1 - b_2.
        return c_0 + scaled_t * b_j - b_i;
      }
  }
}

template<typename Vector>
Vector EvaluationHelper<Vector>::coefficients(int const index) const {
  return coefficients_[index];
}

template<typename Vector>
int EvaluationHelper<Vector>::degree() const {
  return degree_;
}

template<typename Scalar, typename Frame, int rank>
EvaluationHelper<Multivector<Scalar, Frame, rank>>::EvaluationHelper(
    std::vector<Multivector<Scalar, Frame, rank>> const& coefficients,
    int const degree) : degree_(degree) {
  for (auto const& coefficient : coefficients) {
    coefficients_.push_back(coefficient.coordinates() / SIUnit<Scalar>());
  }
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>
EvaluationHelper<Multivector<Scalar, Frame, rank>>::EvaluateImplementation(
    double const scaled_t) const {
  double const two_scaled_t = scaled_t + scaled_t;
  R3Element<double> const c_0 = coefficients_[0];
  switch (degree_) {
    case 0:
      return Multivector<double, Frame, rank>(c_0) * SIUnit<Scalar>();
    case 1:
      return Multivector<double, Frame, rank>(
                 c_0 + scaled_t * coefficients_[1]) * SIUnit<Scalar>();
    default:
      // b_degree   = c_degree.
      R3Element<double> b_i = coefficients_[degree_];
      // b_degree-1 = c_degree-1 + 2 t b_degree.
      R3Element<double> b_j = coefficients_[degree_ - 1] + two_scaled_t * b_i;
      int k = degree_ - 3;
      for (; k >= 1; k -= 2) {
        // b_k+1 = c_k+1 + 2 t b_k+2 - b_k+3.
        R3Element<double> const c_kplus1 = coefficients_[k + 1];
        b_i.x = c_kplus1.x + two_scaled_t * b_j.x - b_i.x;
        b_i.y = c_kplus1.y + two_scaled_t * b_j.y - b_i.y;
        b_i.z = c_kplus1.z + two_scaled_t * b_j.z - b_i.z;
        // b_k   = c_k   + 2 t b_k+1 - b_k+2.
        R3Element<double> const c_k = coefficients_[k];
        b_j.x = c_k.x + two_scaled_t * b_i.x - b_j.x;
        b_j.y = c_k.y + two_scaled_t * b_i.y - b_j.y;
        b_j.z = c_k.z + two_scaled_t * b_i.z - b_j.z;
      }
      if (k == 0) {
        // b_1 = c_1 + 2 t b_2 - b_3.
        b_i = coefficients_[1] + two_scaled_t * b_j - b_i;
        // c_0 + t b_1 - b_2.
        return Multivector<double, Frame, rank>(
                   c_0 + scaled_t * b_i - b_j) * SIUnit<Scalar>();
      } else {
        // c_0 + t b_1 - b_2.
        return Multivector<double, Frame, rank>(
                   c_0 + scaled_t * b_j - b_i) * SIUnit<Scalar>();
      }
    }
}

template<typename Scalar, typename Frame, int rank>
Multivector<Scalar, Frame, rank>
EvaluationHelper<Multivector<Scalar, Frame, rank>>::coefficients(
    int const index) const {
  return Multivector<double, Frame, rank>(
             coefficients_[index]) * SIUnit<Scalar>();
}

template<typename Scalar, typename Frame, int rank>
int EvaluationHelper<Multivector<Scalar, Frame, rank>>::degree() const {
  return degree_;
}

}  // namespace internal

template<typename Vector>
ЧебышёвSeries<Vector>::ЧебышёвSeries(std::vector<Vector> const& coefficients,
                                     Instant const& t_min,
                                     Instant const& t_max)
    : t_min_(t_min),
      t_max_(t_max),
      helper_(coefficients,
              /*degree=*/static_cast<int>(coefficients.size()) - 1) {
  CHECK_LE(0, helper_.degree()) << "Degree must be at least 0";
  CHECK_LT(t_min_, t_max_) << "Time interval must not be empty";
  // Precomputed to save operations at the expense of some accuracy loss.
  Time const duration = t_max_ - t_min_;
  t_mean_ = t_min_ + 0.5 * duration;
  two_over_duration_ = 2 / duration;
}


template<typename Vector>
bool ЧебышёвSeries<Vector>::operator==(ЧебышёвSeries const& right) const {
  if (helper_.degree() != right.helper_.degree()) {
    return false;
  }
  for (int k = 0; k < helper_.degree(); ++k) {
    if (helper_.coefficients(k) != right.helper_.coefficients(k)) {
      return false;
    }
  }
  return t_min_ == right.t_min_ &&
         t_max_ == right.t_max_;
}

template<typename Vector>
bool ЧебышёвSeries<Vector>::operator!=(ЧебышёвSeries const& right) const {
  return !ЧебышёвSeries<Vector>::operator==(right);
}

template<typename Vector>
Instant const& ЧебышёвSeries<Vector>::t_min() const {
  return t_min_;
}

template<typename Vector>
Instant const& ЧебышёвSeries<Vector>::t_max() const {
  return t_max_;
}

template<typename Vector>
Vector ЧебышёвSeries<Vector>::last_coefficient() const {
  return helper_.coefficients(helper_.degree());
}

template<typename Vector>
Vector ЧебышёвSeries<Vector>::Evaluate(Instant const& t) const {
  double const scaled_t = (t - t_mean_) * two_over_duration_;
  // We have to allow |scaled_t| to go slightly out of [-1, 1] because of
  // computation errors.  But if it goes too far, something is broken.
  // TODO(phl): This should use DCHECK but these macros don't work because the
  // Principia projects don't define NDEBUG.
#ifdef _DEBUG
  CHECK_LE(scaled_t, 1.1);
  CHECK_GE(scaled_t, -1.1);
#endif

  return helper_.EvaluateImplementation(scaled_t);
}

template<typename Vector>
Variation<Vector> ЧебышёвSeries<Vector>::EvaluateDerivative(
    Instant const& t) const {
  double const scaled_t = (t - t_mean_) * two_over_duration_;
  double const two_scaled_t = scaled_t + scaled_t;
  // We have to allow |scaled_t| to go slightly out of [-1, 1] because of
  // computation errors.  But if it goes too far, something is broken.
  // TODO(phl): See above.
#ifdef _DEBUG
  CHECK_LE(scaled_t, 1.1);
  CHECK_GE(scaled_t, -1.1);
#endif

  Vector b_kplus2_vector{};
  Vector b_kplus1_vector{};
  Vector* b_kplus2 = &b_kplus2_vector;
  Vector* b_kplus1 = &b_kplus1_vector;
  Vector* const& b_k = b_kplus2;  // An overlay.
  for (int k = helper_.degree() - 1; k >= 1; --k) {
    *b_k = helper_.coefficients(k + 1) * (k + 1) +
           two_scaled_t * *b_kplus1 - *b_kplus2;
    Vector* const last_b_k = b_k;
    b_kplus2 = b_kplus1;
    b_kplus1 = last_b_k;
  }
  return (helper_.coefficients(1) + two_scaled_t * *b_kplus1 - *b_kplus2) *
             two_over_duration_;
}

template<typename Vector>
void ЧебышёвSeries<Vector>::WriteToMessage(
    not_null<serialization::ЧебышёвSeries*> const message) const {
  using Serializer = DoubleOrQuantityOrMultivectorSerializer<
                          Vector,
                          serialization::ЧебышёвSeries::Coefficient>;

  for (int k = 0; k <= helper_.degree(); ++k) {
    Serializer::WriteToMessage(helper_.coefficients(k),
                               message->add_coefficient());
  }
  t_min_.WriteToMessage(message->mutable_t_min());
  t_max_.WriteToMessage(message->mutable_t_max());
}

template<typename Vector>
ЧебышёвSeries<Vector> ЧебышёвSeries<Vector>::ReadFromMessage(
    serialization::ЧебышёвSeries const& message) {
  using Serializer = DoubleOrQuantityOrMultivectorSerializer<
                          Vector,
                          serialization::ЧебышёвSeries::Coefficient>;

  std::vector<Vector> coefficients;
  coefficients.reserve(message.coefficient_size());
  for (auto const& coefficient : message.coefficient()) {
    coefficients.push_back(Serializer::ReadFromMessage(coefficient));
  }
  return ЧебышёвSeries(coefficients,
                       Instant::ReadFromMessage(message.t_min()),
                       Instant::ReadFromMessage(message.t_max()));
}

template<typename Vector>
ЧебышёвSeries<Vector> ЧебышёвSeries<Vector>::NewhallApproximation(
    int const degree,
    std::vector<Vector> const& q,
    std::vector<Variation<Vector>> const& v,
    Instant const& t_min,
    Instant const& t_max) {
  // Only supports 8 divisions for now.
  int const kDivisions = 8;
  CHECK_EQ(kDivisions + 1, q.size());
  CHECK_EQ(kDivisions + 1, v.size());

  Time const duration_over_two = 0.5 * (t_max - t_min);

  // Tricky.  The order in Newhall's matrices is such that the entries for the
  // largest time occur first.
  FixedVector<Vector, 2 * kDivisions + 2> qv;
  for (int i = 0, j = 2 * kDivisions;
       i < kDivisions + 1 && j >= 0;
       ++i, j -= 2) {
    qv[j] = q[i];
    qv[j + 1] = v[i] * duration_over_two;
  }

  std::vector<Vector> coefficients;
  coefficients.reserve(degree);
  switch (degree) {
    case 3:
      coefficients = newhall_c_matrix_degree_3_divisions_8_w04 * qv;
      break;
    case 4:
      coefficients = newhall_c_matrix_degree_4_divisions_8_w04 * qv;
      break;
    case 5:
      coefficients = newhall_c_matrix_degree_5_divisions_8_w04 * qv;
      break;
    case 6:
      coefficients = newhall_c_matrix_degree_6_divisions_8_w04 * qv;
      break;
    case 7:
      coefficients = newhall_c_matrix_degree_7_divisions_8_w04 * qv;
      break;
    case 8:
      coefficients = newhall_c_matrix_degree_8_divisions_8_w04 * qv;
      break;
    case 9:
      coefficients = newhall_c_matrix_degree_9_divisions_8_w04 * qv;
      break;
    case 10:
      coefficients = newhall_c_matrix_degree_10_divisions_8_w04 * qv;
      break;
    case 11:
      coefficients = newhall_c_matrix_degree_11_divisions_8_w04 * qv;
      break;
    case 12:
      coefficients = newhall_c_matrix_degree_12_divisions_8_w04 * qv;
      break;
    case 13:
      coefficients = newhall_c_matrix_degree_13_divisions_8_w04 * qv;
      break;
    case 14:
      coefficients = newhall_c_matrix_degree_14_divisions_8_w04 * qv;
      break;
    case 15:
      coefficients = newhall_c_matrix_degree_15_divisions_8_w04 * qv;
      break;
    case 16:
      coefficients = newhall_c_matrix_degree_16_divisions_8_w04 * qv;
      break;
    case 17:
      coefficients = newhall_c_matrix_degree_17_divisions_8_w04 * qv;
      break;
    default:
      LOG(FATAL) << "Unexpected degree " << degree;
      break;
  }
  CHECK_EQ(degree + 1, coefficients.size());
  return ЧебышёвSeries(coefficients, t_min, t_max);
}

}  // namespace numerics
}  // namespace principia
