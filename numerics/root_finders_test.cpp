
#include <sstream>

#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/point.hpp"
#include "gmock/gmock-generated-matchers.h"
#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "numerics/root_finders.hpp"
#include "physics/discrete_trajectory.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using geometry::Instant;
using quantities::Acceleration;
using quantities::Length;
using quantities::Pow;
using quantities::SIUnit;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Le;

namespace numerics {

class RootFindersTest : public ::testing::Test {};

// Solving Δt * Δt == n.
TEST_F(RootFindersTest, SquareRoots) {
  Instant const t_0;
  Instant const t_max = t_0 + 10 * Second;
  Length const n_max = Pow<2>(t_max - t_0) * SIUnit<Acceleration>();
  for (Length n = 1 * Metre; n < n_max; n += 1 * Metre) {
    int evaluations = 0;
    auto const equation = [t_0, n, &evaluations](Instant const& t) {
      ++evaluations;
      return Pow<2>(t - t_0) * SIUnit<Acceleration>() - n;
    };
    EXPECT_THAT(Bisect(equation, t_0, t_max) - t_0,
                AlmostEquals(Sqrt(n / SIUnit<Acceleration>()), 0, 1));
    if (n == 25 * Metre) {
      EXPECT_EQ(3, evaluations);
    } else {
      EXPECT_THAT(evaluations, AllOf(Ge(49), Le(58)));
    }
  }
}

}  // namespace numerics
}  // namespace principia
