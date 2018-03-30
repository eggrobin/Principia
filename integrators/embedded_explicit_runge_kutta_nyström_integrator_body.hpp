﻿
#pragma once

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <experimental/optional>
#include <vector>

#include "geometry/sign.hpp"
#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace integrators {
namespace internal_embedded_explicit_runge_kutta_nyström_integrator {

using base::make_not_null_unique;
using geometry::Sign;
using numerics::DoublePrecision;
using quantities::DebugString;
using quantities::Difference;
using quantities::Quotient;

template<typename Method, typename Position>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>::
EmbeddedExplicitRungeKuttaNyströmIntegrator()
    : AdaptiveStepSizeIntegrator<
          SpecialSecondOrderDifferentialEquation<Position>>(Method::kind) {
  // the first node is always 0 in an explicit method.
  CHECK_EQ(0.0, c_[0]);
  if (first_same_as_last) {
    // Check that the conditions for the FSAL property are satisfied, see for
    // instance Dormand, El-Mikkawy and Prince (1986),
    // Families of Runge-Kutta-Nyström formulae, equation 3.1.
    CHECK_EQ(1.0, c_[stages_ - 1]);
    CHECK_EQ(0.0, b_hat_[stages_ - 1]);
    for (int j = 0; j < stages_ - 1; ++j) {
      CHECK_EQ(b_hat_[j], a_[stages_ - 1][j]);
    }
  }
}

template<typename Method, typename Position>
Status EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>::
Instance::Solve(Instant const& t_final) {
  using Displacement = typename ODE::Displacement;
  using Velocity = typename ODE::Velocity;
  using Acceleration = typename ODE::Acceleration;

  auto const& a = integrator_.a_;
  auto const& b_hat = integrator_.b_hat_;
  auto const& b_prime_hat = integrator_.b_prime_hat_;
  auto const& b = integrator_.b_;
  auto const& b_prime = integrator_.b_prime_;
  auto const& c = integrator_.c_;

  auto& append_state = this->append_state_;
  auto& current_state = this->current_state_;
  auto& first_use = this->first_use_;
  auto& parameters = this->parameters_;
  auto const& equation = this->equation_;

  // |current_state| gets updated as the integration progresses to allow
  // restartability.

  // State before the last, truncated step.
  std::experimental::optional<typename ODE::SystemState> final_state;

  // Argument checks.
  int const dimension = current_state.positions.size();
  Sign const integration_direction = Sign(parameters.first_time_step);
  if (integration_direction.Positive()) {
    // Integrating forward.
    CHECK_LT(current_state.time.value, t_final);
  } else {
    // Integrating backward.
    CHECK_GT(current_state.time.value, t_final);
  }
  CHECK(first_use || !parameters.last_step_is_exact)
      << "Cannot reuse an instance where the last step is exact";
  first_use = false;

  // Time step.  Updated as the integration progresses to allow restartability.
  Time& h = this->time_step_;
  // Current time.  This is a non-const reference whose purpose is to make the
  // equations more readable.
  DoublePrecision<Instant>& t = current_state.time;

  // Position increment (high-order).
  std::vector<Displacement> Δq_hat(dimension);
  // Velocity increment (high-order).
  std::vector<Velocity> Δv_hat(dimension);
  // Current position.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Position>>& q_hat = current_state.positions;
  // Current velocity.  This is a non-const reference whose purpose is to make
  // the equations more readable.
  std::vector<DoublePrecision<Velocity>>& v_hat = current_state.velocities;

  // Difference between the low- and high-order approximations.
  typename ODE::SystemStateError error_estimate;
  error_estimate.position_error.resize(dimension);
  error_estimate.velocity_error.resize(dimension);

  // Current Runge-Kutta-Nyström stage.
  std::vector<Position> q_stage(dimension);
  // Accelerations at each stage.
  // TODO(egg): this is a rectangular container, use something more appropriate.
  std::vector<std::vector<Acceleration>> g(stages_);
  for (auto& g_stage : g) {
    g_stage.resize(dimension);
  }

  bool at_end = false;
  double tolerance_to_error_ratio;

  // The first stage of the Runge-Kutta-Nyström iteration.  In the FSAL case,
  // |first_stage == 1| after the first step, since the first RHS evaluation has
  // already occured in the previous step.  In the non-FSAL case and in the
  // first step of the FSAL case, |first_stage == 0|.
  int first_stage = 0;

  // The number of steps already performed.
  std::int64_t step_count = 0;

  Status status;

  // No step size control on the first step.  If this instance is being
  // restarted we already have a value of |h| suitable for the next step, based
  // on the computation of |tolerance_to_error_ratio_| during the last
  // invocation.
  goto runge_kutta_nyström_step;

  while (!at_end) {
    // Compute the next step with decreasing step sizes until the error is
    // tolerable.
    do {
      // Adapt step size.
      // TODO(egg): find out whether there's a smarter way to compute that root,
      // especially since we make the order compile-time.
      h *= parameters.safety_factor *
               std::pow(tolerance_to_error_ratio, 1.0 / (lower_order + 1));
      // TODO(egg): should we check whether it vanishes in double precision
      // instead?
      if (t.value + (t.error + h) == t.value) {
        return Status(termination_condition::VanishingStepSize,
                      "At time " + DebugString(t.value) +
                          ", step size is effectively zero.  Singularity or "
                          "stiff system suspected.");
      }

    runge_kutta_nyström_step:
      // Termination condition.
      if (parameters.last_step_is_exact) {
        Time const time_to_end = (t_final - t.value) - t.error;
        at_end = integration_direction * h >=
                 integration_direction * time_to_end;
        if (at_end) {
          // The chosen step size will overshoot.  Clip it to just reach the
          // end, and terminate if the step is accepted.
          h = time_to_end;
          final_state = current_state;
        }
      }

      // Runge-Kutta-Nyström iteration; fills |g|.
      for (int i = first_stage; i < stages_; ++i) {
        Instant const t_stage = t.value + (t.error + c[i] * h);
        for (int k = 0; k < dimension; ++k) {
          Acceleration Σj_a_ij_g_jk{};
          for (int j = 0; j < i; ++j) {
            Σj_a_ij_g_jk += a[i][j] * g[j][k];
          }
          q_stage[k] = q_hat[k].value +
                           h * (c[i] * v_hat[k].value + h * Σj_a_ij_g_jk);
        }
        status.Update(equation.compute_acceleration(t_stage, q_stage, g[i]));
      }

      // Increment computation and step size control.
      for (int k = 0; k < dimension; ++k) {
        Acceleration Σi_b_hat_i_g_ik{};
        Acceleration Σi_b_i_g_ik{};
        Acceleration Σi_b_prime_hat_i_g_ik{};
        Acceleration Σi_b_prime_i_g_ik{};
        // Please keep the eight assigments below aligned, they become illegible
        // otherwise.
        for (int i = 0; i < stages_; ++i) {
          Σi_b_hat_i_g_ik       += b_hat[i] * g[i][k];
          Σi_b_i_g_ik           += b[i] * g[i][k];
          Σi_b_prime_hat_i_g_ik += b_prime_hat[i] * g[i][k];
          Σi_b_prime_i_g_ik     += b_prime[i] * g[i][k];
        }
        // The hat-less Δq and Δv are the low-order increments.
        Δq_hat[k]               = h * (h * (Σi_b_hat_i_g_ik) + v_hat[k].value);
        Displacement const Δq_k = h * (h * (Σi_b_i_g_ik) + v_hat[k].value);
        Δv_hat[k]               = h * Σi_b_prime_hat_i_g_ik;
        Velocity const Δv_k     = h * Σi_b_prime_i_g_ik;

        error_estimate.position_error[k] = Δq_k - Δq_hat[k];
        error_estimate.velocity_error[k] = Δv_k - Δv_hat[k];
      }
      tolerance_to_error_ratio =
          this->tolerance_to_error_ratio_(h, error_estimate);
    } while (tolerance_to_error_ratio < 1.0);

    if (!parameters.last_step_is_exact && t.value + (t.error + h) > t_final) {
      // We did overshoot.  Drop the point that we just computed and exit.
      final_state = current_state;
      break;
    }

    if (first_same_as_last) {
      using std::swap;
      swap(g.front(), g.back());
      first_stage = 1;
    }

    // Increment the solution with the high-order approximation.
    t.Increment(h);
    for (int k = 0; k < dimension; ++k) {
      q_hat[k].Increment(Δq_hat[k]);
      v_hat[k].Increment(Δv_hat[k]);
    }
    append_state(current_state);
    ++step_count;
    if (step_count == parameters.max_steps && !at_end) {
      return Status(termination_condition::ReachedMaximalStepCount,
                    "Reached maximum step count " +
                        std::to_string(parameters.max_steps) +
                        " at time " + DebugString(t.value) +
                        "; requested t_final is " + DebugString(t_final) +
                        ".");
    }
  }
  // The resolution is restartable from the last non-truncated state.
  CHECK(final_state);
  current_state = *final_state;
  return status;
}

template<typename Method, typename Position>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position> const&
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>::
Instance::integrator() const {
  return integrator_;
}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    SpecialSecondOrderDifferentialEquation<Position>>::Instance>>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>::
Instance::Clone() const {
  return std::unique_ptr<Instance>(new Instance(*this));
}

template<typename Method, typename Position>
void EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>::
Instance::WriteToMessage(
    not_null<serialization::IntegratorInstance*> message) const {
  AdaptiveStepSizeIntegrator<ODE>::Instance::WriteToMessage(message);
  auto* const extension =
      message
          ->MutableExtension(
              serialization::AdaptiveStepSizeIntegratorInstance::extension)
          ->MutableExtension(
              serialization::
                  EmbeddedExplicitRungeKuttaNystromIntegratorInstance::
                      extension);
}

template<typename Method, typename Position>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>::
Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio,
    Parameters const& parameters,
    EmbeddedExplicitRungeKuttaNyströmIntegrator const& integrator)
    : AdaptiveStepSizeIntegrator<ODE>::Instance(
          problem, append_state, tolerance_to_error_ratio, parameters),
      integrator_(integrator) {}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    SpecialSecondOrderDifferentialEquation<Position>>::Instance>>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>::
NewInstance(IntegrationProblem<ODE> const& problem,
            AppendState const& append_state,
            ToleranceToErrorRatio const& tolerance_to_error_ratio,
            Parameters const& parameters) const {
  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<Instance>(new Instance(problem,
                                                append_state,
                                                tolerance_to_error_ratio,
                                                parameters,
                                                *this));
}

template<typename Method, typename Position>
not_null<std::unique_ptr<typename Integrator<
    SpecialSecondOrderDifferentialEquation<Position>>::Instance>>
EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position>::
ReadFromMessage(
    serialization::AdaptiveStepSizeIntegratorInstance const& message,
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    ToleranceToErrorRatio const& tolerance_to_error_ratio,
    Parameters const& parameters) const {
  CHECK(message.HasExtension(
      serialization::EmbeddedExplicitRungeKuttaNystromIntegratorInstance::
          extension))
      << message.DebugString();

  // Cannot use |make_not_null_unique| because the constructor of |Instance| is
  // private.
  return std::unique_ptr<typename AdaptiveStepSizeIntegrator<ODE>::Instance>(
      new Instance(problem,
                   append_state,
                   tolerance_to_error_ratio,
                   parameters,
                   *this));
}

}  // namespace internal_embedded_explicit_runge_kutta_nyström_integrator

template<typename Method, typename Position>
internal_embedded_explicit_runge_kutta_nyström_integrator::
    EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position> const&
EmbeddedExplicitRungeKuttaNyströmIntegrator() {
  static_assert(
      std::is_base_of<methods::EmbeddedExplicitRungeKuttaNyström,
                      Method>::value,
      "Method must be derived from EmbeddedExplicitRungeKuttaNyström");
  static internal_embedded_explicit_runge_kutta_nyström_integrator::
      EmbeddedExplicitRungeKuttaNyströmIntegrator<Method, Position> const
          integrator;
  return integrator;
}

}  // namespace integrators
}  // namespace principia
