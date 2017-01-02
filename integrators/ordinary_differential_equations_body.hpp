﻿
#pragma once

#include "integrators/ordinary_differential_equations.hpp"

#include <vector>

#include "base/macros.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"

namespace principia {
namespace integrators {
namespace internal_ordinary_differential_equations {

// TODO(egg): for some mysterious reason MSVC wants the full
// |typename SpecialSecondOrderDifferentialEquation<Position_>::Position|
// where |Position| would be enough.
template<typename Position_>
SpecialSecondOrderDifferentialEquation<Position_>::SystemState::SystemState(
    std::vector<typename SpecialSecondOrderDifferentialEquation<
        Position_>::Position> const& q,
    std::vector<typename SpecialSecondOrderDifferentialEquation<
        Position_>::Velocity> const& v,
    Instant const& t)
    : time(t) {
  for (int i = 0; i < q.size(); ++i) {
    positions.emplace_back(q[i]);
    velocities.emplace_back(v[i]);
  }
}

template<typename Position>
void
SpecialSecondOrderDifferentialEquation<Position>::SystemState::WriteToMessage(
        not_null<serialization::SystemState*> const message) const {
  for (auto const& position : positions) {
    position.WriteToMessage(message->add_position());
  }
  for (auto const& velocity : velocities) {
    velocity.WriteToMessage(message->add_velocity());
  }
  time.WriteToMessage(message->mutable_time());
}

template<typename Position>
typename SpecialSecondOrderDifferentialEquation<Position>::SystemState
SpecialSecondOrderDifferentialEquation<Position>::SystemState::ReadFromMessage(
        serialization::SystemState const& message) {
  SystemState system_state;
  for (auto const p : message.position()) {
    system_state.positions.push_back(
        DoublePrecision<Position>::ReadFromMessage(p));
  }
  for (auto const v : message.velocity()) {
    system_state.velocities.push_back(
        DoublePrecision<Velocity>::ReadFromMessage(v));
  }
  system_state.time = DoublePrecision<Instant>::ReadFromMessage(message.time());
  return system_state;
}

template<typename DifferentialEquation>
Integrator<DifferentialEquation>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state)
    : equation_(problem.equation),
      current_state_(*problem.initial_state),
      append_state_(std::move(append_state)) {
  CHECK_EQ(current_state_.positions.size(),
           current_state_.velocities.size());
}

template<typename DifferentialEquation>
Instant const& Integrator<DifferentialEquation>::Instance::time() const {
  return current_state_.time.value;
}

template<typename DifferentialEquation>
Integrator<DifferentialEquation>::Instance::Instance() {}

template<typename DifferentialEquation>
FixedStepSizeIntegrator<DifferentialEquation>::FixedStepSizeIntegrator(
    serialization::FixedStepSizeIntegrator::Kind const kind) : kind_(kind) {}

template<typename DifferentialEquation>
void FixedStepSizeIntegrator<DifferentialEquation>::WriteToMessage(
    not_null<serialization::FixedStepSizeIntegrator*> const message) const {
  message->set_kind(kind_);
}

template<typename DifferentialEquation>
FixedStepSizeIntegrator<DifferentialEquation>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    Time const& step)
    : Integrator<ODE>::Instance(problem, std::move(append_state)),
      step_(step) {
  CHECK_NE(Time(), step_);
}

template<typename DifferentialEquation>
FixedStepSizeIntegrator<DifferentialEquation> const&
FixedStepSizeIntegrator<DifferentialEquation>::ReadFromMessage(
      serialization::FixedStepSizeIntegrator const& message) {
  using FSSI = serialization::FixedStepSizeIntegrator;
  switch (message.kind()) {
    case FSSI::BLANES_MOAN_2002_SRKN_6B:
      return BlanesMoan2002SRKN6B<typename DifferentialEquation::Position>();
    case FSSI::BLANES_MOAN_2002_SRKN_11B:
      return BlanesMoan2002SRKN11B<typename DifferentialEquation::Position>();
    case FSSI::BLANES_MOAN_2002_SRKN_14A:
      return BlanesMoan2002SRKN14A<typename DifferentialEquation::Position>();
    case FSSI::MCLACHLAN_1995_SB3A_4:
      return McLachlan1995SB3A4<typename DifferentialEquation::Position>();
    case FSSI::MCLACHLAN_1995_SB3A_5:
      return McLachlan1995SB3A5<typename DifferentialEquation::Position>();
    case FSSI::MCLACHLAN_ATELA_1992_ORDER_4_OPTIMAL:
      return McLachlanAtela1992Order4Optimal<
                 typename DifferentialEquation::Position>();
    case FSSI::MCLACHLAN_ATELA_1992_ORDER_5_OPTIMAL:
      return McLachlanAtela1992Order5Optimal<
                 typename DifferentialEquation::Position>();
    case FSSI::OKUNBOR_SKEEL_1994_ORDER_6_METHOD_13:
      return OkunborSkeel1994Order6Method13<
                 typename DifferentialEquation::Position>();
    default:
      LOG(FATAL) << message.kind();
      base::noreturn();
  }
}

template<typename DifferentialEquation>
AdaptiveStepSizeIntegrator<DifferentialEquation>::AdaptiveStepSizeIntegrator(
    serialization::AdaptiveStepSizeIntegrator::Kind const kind) : kind_(kind) {}

template <typename DifferentialEquation>
AdaptiveStepSizeIntegrator<DifferentialEquation>::Instance::Instance(
    IntegrationProblem<ODE> const& problem,
    AppendState const& append_state,
    AdaptiveStepSize<ODE> const& adaptive_step_size)
    : Integrator<ODE>::Instance(problem, std::move(append_state)),
      adaptive_step_size_(adaptive_step_size) {
  CHECK_NE(Time(), adaptive_step_size.first_time_step);
  CHECK_GT(adaptive_step_size.safety_factor, 0);
  CHECK_LT(adaptive_step_size.safety_factor, 1);
}

template<typename DifferentialEquation>
void AdaptiveStepSizeIntegrator<DifferentialEquation>::WriteToMessage(
    not_null<serialization::AdaptiveStepSizeIntegrator*> const message) const {
  message->set_kind(kind_);
}

template<typename DifferentialEquation>
AdaptiveStepSizeIntegrator<DifferentialEquation> const&
AdaptiveStepSizeIntegrator<DifferentialEquation>::ReadFromMessage(
    serialization::AdaptiveStepSizeIntegrator const& message) {
  using ASSI = serialization::AdaptiveStepSizeIntegrator;
  switch (message.kind()) {
    case ASSI::DORMAND_ELMIKKAWY_PRINCE_1986_RKN_434FM:
      return DormandElMikkawyPrince1986RKN434FM<
                 typename DifferentialEquation::Position>();
    default:
      LOG(FATAL) << message.kind();
      base::noreturn();
  }
}

}  // namespace internal_ordinary_differential_equations
}  // namespace integrators
}  // namespace principia
