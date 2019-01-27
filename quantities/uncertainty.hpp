#pragma once

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {

// In this file:
// — GUM refers to:
//   Guide to the expression of uncertainty in measurement, first edition,
//   JCGM 100:2008;
// — VIM refers to:
//   Vocabulaire international de métrologie, 3ᵉ édition, JCGM 200:2012.

// A measurement result, expressed as a measured quantity value and a
// measurement uncertainty (VIM 2.9 note 2).
// The measurement uncertainty is expressed as a standard deviation, i.e., it is
// a standard uncertainty (VIM 2.30, GUM 2.3.1).
// T may be a |Quantity| or a |Point<Quantity>|.
template<typename T>
struct MeasurementResult {
  T measured_value;
  Difference<T> standard_uncertainty;
};

std::ostream& operator<<(std::ostream& out,
                         MeasurementResult<double> measurement);

// Multiplication or division of a measurement result by a an exact value.
template<typename T, typename U>
MeasurementResult<Product<T, U>> operator*(MeasurementResult<T> const& left,
                                           U const& right);
template<typename T, typename U>
MeasurementResult<Product<T, U>> operator*(T const& left,
                                           MeasurementResult<U> const& right);
template<typename T, typename U>
MeasurementResult<Quotient<T, U>> operator/(MeasurementResult<T> const& left,
                                            U const& right);

// Inverse of a measurement result; the uncertainty is propagated according to
// GUM 5.1.2.
template<typename T, typename U>
MeasurementResult<Quotient<T, U>> operator/(T const& left,
                                            MeasurementResult<U> const& right);

// The following functions provide type A evaluations of standard uncertainty
// (VIM 2.28, GUM 4.2), for various measurement models (VIM 2.48).

// The output quantities are obtained by projections of the input quantities,
// which we express in the formalism of Grant Foster (1996), Time series
// analysis by projection.

// A quantity estimated as an average of independent observations.
// |measured_values| must be independent and identically distributed.
// In the formalism of Foster (1996), this is the projection onto |𝟏⟩ with the
// uncertainty evaluated under the basic null hypothesis.
template<typename T>
MeasurementResult<T> AverageOfIndependent(
    std::vector<T> const& measured_values);

// An estimation of the mean of an indexed time series.
// The values of the time series must be identically distributed, but need not
// be independent.  In practice, this allows for the computation of the mean in
// the presence of a strong signal without having to model the signal.
// In the formalism of Foster (1996), this is the projection onto |𝟏⟩ with the
// uncertainty evaluated under the correlated random variable hypothesis.
template<typename T>
MeasurementResult<T> AverageOfCorrelated(std::vector<T> const& time_series);

// The result of a linear regression; in the formalism of Foster (1996), this is
// the projection onto the basis consisting of:
// — |𝟏⟩, the constant function, and
// — |t - ⟨𝟏|t⟩|𝟏⟩⟩, the zero-mean linear function of the argument t,
// where the uncertainties are estimated under the strong signal hypothesis.
// The uncertainties of |mean_value| and |slope| are independent (this is
// because the basis is orthogonal; see also GUM H.3.5).
template<typename Argument, typename Value>
struct LinearModel {
  // The coefficient of |𝟏⟩.
  MeasurementResult<Value> mean_value;
  // The coefficient of |t - ⟨𝟏|t⟩|𝟏⟩⟩.
  MeasurementResult<Derivative<Value, Argument>> slope;
};

template<typename Argument, typename Value>
LinearModel<Argument, Value> LinearRegression(
    std::vector<Argument> const& arguments,
    std::vector<Value> const& values);

}  // namespace quantities
}  // namespace principia

#include "quantities/uncertainty_body.hpp"
