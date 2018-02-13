
#include <intrin.h>
#include <vector>

#include "benchmark/benchmark.h"

struct V {
  __m128d xy;
  __m128d zt;

  double hadd_norm2() {
    __m128d const zero = _mm_setzero_pd();
    __m128d const z_0 = _mm_unpacklo_pd(zt, zero);
    __m128d const z²_0 = _mm_mul_sd(z_0, z_0);
    __m128d const x²_y² = _mm_mul_pd(xy, xy);
    __m128d const x²z²_y² = _mm_add_sd(x²_y², z²_0);
    __m128d const x²y²z²_0 = _mm_hadd_pd(x²z²_y², zero);
    return x²y²z²_0.m128d_f64[0];
  }

  double nohadd_norm2() {
    __m128d const z²_0 = _mm_mul_sd(zt, zt);
    __m128d const x²_y² = _mm_mul_pd(xy, xy);
    __m128d const y²_y² = _mm_castps_pd(
        _mm_movehl_ps(_mm_castpd_ps(x²_y²), _mm_castpd_ps(x²_y²)));
    __m128d const x²y²_y² = _mm_add_sd(x²_y², y²_y²);
    __m128d const x²y²z²_y² = _mm_add_sd(x²y²_y², z²_0);
    return x²y²z²_y².m128d_f64[0];
  }

  double norm2() {
    return xy.m128d_f64[0] * xy.m128d_f64[0] +
           xy.m128d_f64[1] * xy.m128d_f64[1] +
           zt.m128d_f64[0] * zt.m128d_f64[0];
  }
};

#define TENFOLD(x) x; x; x; x; x; x; x; x; x; x

#define DO_TEN(tens, f) \
    results[tens##0] = vs[tens##0].f(); \
    results[tens##1] = vs[tens##1].f(); \
    results[tens##2] = vs[tens##2].f(); \
    results[tens##3] = vs[tens##3].f(); \
    results[tens##4] = vs[tens##4].f(); \
    results[tens##5] = vs[tens##5].f(); \
    results[tens##6] = vs[tens##6].f(); \
    results[tens##7] = vs[tens##7].f(); \
    results[tens##8] = vs[tens##8].f(); \
    results[tens##9] = vs[tens##9].f();

#define BM(f)                                                     \
  \
void BM_##f##_indep(benchmark::State& state) {                    \
    V v = {{{1, 1}}, {{1, 0}}};                                   \
    std::vector<V> vs;                                            \
    std::vector<double> results;                                  \
    for (int i = 0; i < 100; ++i) {                               \
      vs.push_back(v);                                            \
      results.push_back(0);                                       \
    }                                                             \
    for (auto _ : state) {                                        \
      int i = 0;                                                  \
      TENFOLD(TENFOLD(results[i] = vs[i].f(); ++i));              \
    }                                                             \
    state.SetLabel(std::to_string(results[42]));                  \
  \
}                                                            \
  \
BENCHMARK(BM_##f##_indep);                                        \
  \
void BM_##f##_seq(benchmark::State& state) {                      \
    V v = {{{1, 1}}, {{1, 0}}};                                   \
    std::vector<V> vs;                                            \
    double result = 1;                                            \
    for (int i = 0; i < 100; ++i) {                               \
      vs.push_back(v);                                            \
    }                                                             \
    for (auto _ : state) {                                        \
      int i = 0;                                                  \
      TENFOLD(TENFOLD(vs[i].xy = _mm_load1_pd(&result);           \
                      vs[i].zt = _mm_move_sd(vs[i].zt, vs[i].xy); \
                      result = vs[i].f();                         \
                      ++i));                                      \
    }                                                             \
    state.SetLabel(std::to_string(result));                       \
  \
}                                                            \
  \
BENCHMARK(BM_##f##_seq);

BM(norm2);
BM(nohadd_norm2);
BM(hadd_norm2);


