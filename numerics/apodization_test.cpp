
#include "numerics/apodization.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace numerics {

using geometry::Instant;
using quantities::Sqrt;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::VanishesBefore;

class ApodizationTest : public ::testing::Test {
 protected:
  ApodizationTest()
      : t1_(t0_ - 1 * Second),
        t2_(t0_ + 2 * Second),
        mid_(t0_ + 0.5 * Second) {}

  Instant const t0_;
  Instant const t1_;
  Instant const t2_;
  Instant const mid_;
};

TEST_F(ApodizationTest, Dirichlet) {
  auto a = apodization::Dirichlet<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), AlmostEquals(1, 0));
  EXPECT_THAT(a(t0_), AlmostEquals(1, 0));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a(t2_), AlmostEquals(1, 0));
}

TEST_F(ApodizationTest, Sine) {
  auto a = apodization::Sine<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), VanishesBefore(1, 0));
  EXPECT_THAT(a(t0_), AlmostEquals(Sqrt(3) / 2, 1));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a(t2_), VanishesBefore(1, 0));
}

TEST_F(ApodizationTest, Hann) {
  auto a = apodization::Hann<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), AlmostEquals(0, 0));
  EXPECT_THAT(a(t0_), AlmostEquals(0.75, 0));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a(t2_), VanishesBefore(1, 0));
}

TEST_F(ApodizationTest, Hamming) {
  auto a = apodization::Hamming<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), AlmostEquals(2.0 / 23.0, 0));
  EXPECT_THAT(a(t0_), AlmostEquals(71.0 / 92.0, 0));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a(t2_), AlmostEquals(2.0 / 23.0, 0));
}

TEST_F(ApodizationTest, Blackman) {
  auto a = apodization::Blackman<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), VanishesBefore(1, 0));
  EXPECT_THAT(a(t0_), AlmostEquals(0.63, 1));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 1));
  EXPECT_THAT(a(t2_), VanishesBefore(1, 0));
}

TEST_F(ApodizationTest, ExactBlackman) {
  auto a = apodization::ExactBlackman<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), AlmostEquals(8.0 / 1163.0, 37));
  EXPECT_THAT(a(t0_), AlmostEquals(11843.0 / 18608.0, 2));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a(t2_), AlmostEquals(8.0 / 1163.0, 37));
}

TEST_F(ApodizationTest, Nuttall) {
  auto a = apodization::Nuttall<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), VanishesBefore(1, 0));
  EXPECT_THAT(a(t0_), AlmostEquals(0.514746, 1));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a(t2_), VanishesBefore(1, 0));
}

TEST_F(ApodizationTest, BlackmanNuttall) {
  auto a = apodization::BlackmanNuttall<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), AlmostEquals(0.0003628, 703));
  EXPECT_THAT(a(t0_), AlmostEquals(0.5292298, 1));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a(t2_), AlmostEquals(0.0003628, 703));
}

TEST_F(ApodizationTest, BlackmanHarris) {
  auto a = apodization::BlackmanHarris<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), AlmostEquals(0.00006, 151));
  EXPECT_THAT(a(t0_), AlmostEquals(0.520575, 1));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 0));
  EXPECT_THAT(a(t2_), AlmostEquals(0.00006, 151));
}

TEST_F(ApodizationTest, ISO18431_2) {
  auto a = apodization::ISO18431_2<HornerEvaluator>(t1_, t2_);
  EXPECT_THAT(a(t1_), AlmostEquals(-0.00195313 / 4.63867187, 440));
  EXPECT_THAT(a(t0_), AlmostEquals(0.9194336 / 4.63867187, 5));
  EXPECT_THAT(a(mid_), AlmostEquals(1, 1));
  EXPECT_THAT(a(t2_), AlmostEquals(-0.00195313 / 4.63867187, 440));
}

}  // namespace numerics
}  // namespace principia
