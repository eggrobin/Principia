﻿
#include <cmath>
#include <ctime>
#include <memory>
#include <vector>

// NOTE(phl): The glog operations may not work with Unity/Mono.


namespace principia {
namespace integrators {

namespace {

template <typename T>
inline std::vector<T>* PointerOrNew(int const dimension,
                                    std::vector<T>* const in) {
  if (in == nullptr) {
    return new std::vector<T>(dimension);
  } else {
    return in;
  }
}

}  // namespace

inline std::vector<std::vector<double>> const&
SPRKIntegrator::Order5Optimal() const {
  static std::vector<std::vector<double>> const order_5_optimal = {
      { 0.339839625839110000,
       -0.088601336903027329,
        0.5858564768259621188,
       -0.603039356536491888,
        0.3235807965546976394,
        0.4423637942197494587},
      { 0.1193900292875672758,
        0.6989273703824752308,
       -0.1713123582716007754,
        0.4012695022513534480,
        0.0107050818482359840,
       -0.0589796254980311632}};
  return order_5_optimal;
}

inline void SPRKIntegrator::Solve(
      RightHandSideComputation const compute_force,
      AutonomousRightHandSideComputation const compute_velocity,
      Parameters const& parameters,
      Solution* solution) {

  std::vector<double> const& a = parameters.coefficients[0];
  std::vector<double> const& b = parameters.coefficients[1];
  int const stages = b.size();
  int const dimension = parameters.q0.size();

  // Runge-Kutta time weights.
  std::vector<double> c(stages);
  c[0] = 0.0;
  for (int j = 1; j < stages; ++j) {
    c[j] = c[j - 1] + b[j - 1];
  }

  std::unique_ptr<std::vector<double>> q_error(
      PointerOrNew(dimension, parameters.q_error));
  std::unique_ptr<std::vector<double>> p_error(
      PointerOrNew(dimension, parameters.p_error));
  double t_error = parameters.t_error;

  std::vector<std::vector<double>> Δqstages(stages + 1);
  std::vector<std::vector<double>> Δpstages(stages + 1);

  for (int i = 0; i < stages + 1; ++i) {
    Δqstages[i].resize(dimension);
    Δpstages[i].resize(dimension);
  }

  // Result goes here.
  int const capacity = parameters.sampling_period == 0 ?
    1 :
    static_cast<int>(
        ceil((((parameters.tmax - parameters.t0) / parameters.Δt) + 1) /
                parameters.sampling_period)) + 1;
  std::vector<std::vector<double>> q;
  q.reserve(capacity);
  std::vector<std::vector<double>> p;
  p.reserve(capacity);
  std::vector<double> t;
  t.reserve(capacity);
  if (parameters.sampling_period != 0) {
    q.push_back(parameters.q0);
    p.push_back(parameters.p0);
    t.push_back(parameters.t0);
  }

  std::vector<double> q_last(parameters.q0);
  std::vector<double> p_last(parameters.p0);
  double t_last = parameters.t0;
  int sampling_phase = 0;

  std::vector<double> q_stage(dimension);
  std::vector<double> p_stage(dimension);
  double tn = parameters.t0;  // Current time.
  double const h = parameters.Δt;  // Constant for now.
  std::vector<double> f(dimension);  // Current forces.
  std::vector<double> v(dimension);  // Current velocities.

#ifdef TRACE
  int percentage = 0;
  clock_t running_time = clock();
#endif

  // Integration.  For details see Wolfram Reference,
  // http://reference.wolfram.com/mathematica/tutorial/NDSolveSPRK.html#74387056
  while (tn < parameters.tmax) {
    // Increment SPRK step from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 3.
    for (int k = 0; k < dimension; ++k) {
      Δqstages[0][k] = 0;
      Δpstages[0][k] = 0;
      q_stage[k] = q_last[k];
    }
    for (int i = 1; i < stages + 1; ++i) {
      // Beware, the p/q order matters here, the two computations depend on one
      // another.
      compute_force(tn + c[i - 1] * h, q_stage, &f);
      for (int k = 0; k < dimension; ++k) {
        Δpstages[i][k] = Δpstages[i - 1][k] + h * b[i - 1] * f[k];
        p_stage[k] = p_last[k] + Δpstages[i][k];
      }
      compute_velocity(p_stage, &v);
      for (int k = 0; k < dimension; ++k) {
        Δqstages[i][k] = Δqstages[i - 1][k] + h * a[i - 1] * v[k];
        q_stage[k] = q_last[k] + Δqstages[i][k];
      }
    }
    // Compensated summation from "'SymplecticPartitionedRungeKutta' Method
    // for NDSolve", algorithm 2.
    for (int k = 0; k < dimension; ++k) {
      double const Δq = Δqstages[stages][k] + (*q_error)[k];
      q_stage[k] = q_last[k] + Δq;
      (*q_error)[k] = (q_last[k] - q_stage[k]) + Δq;
      q_last[k] = q_stage[k];
      double const Δp = Δpstages[stages][k] + (*p_error)[k];
      p_stage[k] = p_last[k] + Δp;
      (*p_error)[k] = (p_last[k] - p_stage[k]) + Δp;
      p_last[k] = p_stage[k];
    }

    double const δt = h + t_error;
    tn += δt;
    t_error = (t_last - tn) + δt;
    t_last = tn;

    if (parameters.sampling_period != 0) {
      if (sampling_phase % parameters.sampling_period == 0) {
        t.push_back(tn);
        q.push_back(q_stage);
        p.push_back(p_stage);
      }
      ++sampling_phase;
    }

#ifdef TRACE
    running_time += clock();
    if (floor(tn / parameters.tmax * 100) > percentage) {
      LOG(ERROR) << "SPRK: " << percentage << "%\ttn = " << tn
                 <<"\tRunning time: " << running_time / (CLOCKS_PER_SEC / 1000)
                 << " ms";
      ++percentage;
    }
    running_time -= clock();
#endif
  }
  if (parameters.sampling_period == 0) {
    t.push_back(tn);
    q.push_back(q_stage);
    p.push_back(p_stage);
  }

  solution->momentum.resize(dimension);
  solution->position.resize(dimension);
  for (size_t i = 0; i < p.size(); ++i) {
    for (int j = 0; j < dimension; ++j) {
      if (i == 0) {
        solution->position[j].quantities.reserve(q.size());
        solution->momentum[j].quantities.reserve(p.size());
      }
      solution->position[j].quantities.push_back(q[i][j]);
      solution->momentum[j].quantities.push_back(p[i][j]);
    }
  }
  for (int j = 0; j < dimension; ++j) {
    solution->position[j].error = (*q_error)[j];
    solution->momentum[j].error = (*p_error)[j];
  }
  solution->time.quantities = t;
  solution->time.error = t_error;
}

}  // namespace integrators
}  // namespace principia

#undef TRACE
