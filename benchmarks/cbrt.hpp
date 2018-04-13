#pragma once

#include <emmintrin.h>

#include <map>
#include <string>

namespace principia {
namespace numerics {

using SingleParameterScalarFunction = double(double);

class CubeRootRegistry {
 public:
  static CubeRootRegistry& Instance();

  SingleParameterScalarFunction* Register(std::string const& name,
                                          SingleParameterScalarFunction* cbrt);

  std::map<std::string, SingleParameterScalarFunction*> const& methods() const;

 private:
  CubeRootRegistry() = default;
  std::map<std::string, double (*)(double)> methods_;
};

inline std::uint64_t to_integer(double x) {
  return _mm_cvtsi128_si64(_mm_castpd_si128(_mm_set_sd(x)));
}

inline double to_double(std::uint64_t x) {
  return _mm_cvtsd_f64(_mm_castsi128_pd(_mm_cvtsi64_si128(x)));
}

struct RoundedReal {
  double nearest_rounding;
  double furthest_rounding;
  double nearest_ulps;
};

RoundedReal correct_cube_root(double const y);

#define PRINCIPIA_REGISTER_CBRT(name)                                 \
  static principia::numerics::SingleParameterScalarFunction* const    \
      name##_cbrt =                                                   \
          principia::numerics::CubeRootRegistry::Instance().Register( \
              #name, &name::cbrt)

}  // namespace numerics
}  // namespace principia
