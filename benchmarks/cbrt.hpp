#pragma once

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

#define PRINCIPIA_REGISTER_CBRT(name)                                 \
  static principia::numerics::SingleParameterScalarFunction* const    \
      name##_cbrt =                                                   \
          principia::numerics::CubeRootRegistry::Instance().Register( \
              #name, &name::cbrt)

}  // namespace numerics
}  // namespace principia
