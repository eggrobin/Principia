﻿
#include "testing_utilities/approximate_quantity.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace testing_utilities {

using quantities::Area;
using quantities::Frequency;
using quantities::Length;
using quantities::Speed;
using quantities::si::Metre;
using quantities::si::Second;

TEST(ApproximateQuantityTest, Literals_⑴) {
  ApproximateQuantity<double> const l1 = 123.45_⑴;
  EXPECT_THAT(l1.min(), AlmostEquals(123.44, 0));
  EXPECT_THAT(l1.max(), AlmostEquals(123.46, 1));

  ApproximateQuantity<double> const l2 = 123.45e-3_⑴;
  EXPECT_THAT(l2.min(), AlmostEquals(123.44e-3, 1));
  EXPECT_THAT(l2.max(), AlmostEquals(123.46e-3, 0));

  ApproximateQuantity<double> const l3 = 123e3_⑴;
  EXPECT_THAT(l3.min(), AlmostEquals(122e3, 0));
  EXPECT_THAT(l3.max(), AlmostEquals(124e3, 0));

  ApproximateQuantity<double> const l4 = 0x1E3.45p0_⑴;
  EXPECT_THAT(l4.min(), AlmostEquals(0x1E3.44p0, 0));
  EXPECT_THAT(l4.max(), AlmostEquals(0x1E3.46p0, 0));

  ApproximateQuantity<double> const l5 = 123_⑴;
  EXPECT_THAT(l5.min(), AlmostEquals(122, 0));
  EXPECT_THAT(l5.max(), AlmostEquals(124, 0));
}

TEST(ApproximateQuantityTest, Literals_⑵_⑼) {
  ApproximateQuantity<double> const l2 = 123.45_⑵;
  EXPECT_THAT(l2.min(), AlmostEquals(123.43, 0));
  EXPECT_THAT(l2.max(), AlmostEquals(123.47, 0));

  ApproximateQuantity<double> const l3 = 123.45_⑶;
  EXPECT_THAT(l3.min(), AlmostEquals(123.42, 0));
  EXPECT_THAT(l3.max(), AlmostEquals(123.48, 0));

  ApproximateQuantity<double> const l4 = 123.45_⑷;
  EXPECT_THAT(l4.min(), AlmostEquals(123.41, 0));
  EXPECT_THAT(l4.max(), AlmostEquals(123.49, 1));

  ApproximateQuantity<double> const l5 = 123.45_⑸;
  EXPECT_THAT(l5.min(), AlmostEquals(123.40, 0));
  EXPECT_THAT(l5.max(), AlmostEquals(123.50, 0));

  ApproximateQuantity<double> const l6 = 123.45_⑹;
  EXPECT_THAT(l6.min(), AlmostEquals(123.39, 0));
  EXPECT_THAT(l6.max(), AlmostEquals(123.51, 0));

  ApproximateQuantity<double> const l7 = 123.45_⑺;
  EXPECT_THAT(l7.min(), AlmostEquals(123.38, 1));
  EXPECT_THAT(l7.max(), AlmostEquals(123.52, 0));

  ApproximateQuantity<double> const l8 = 123.45_⑻;
  EXPECT_THAT(l8.min(), AlmostEquals(123.37, 0));
  EXPECT_THAT(l8.max(), AlmostEquals(123.53, 0));

  ApproximateQuantity<double> const l9 = 123.45_⑼;
  EXPECT_THAT(l9.min(), AlmostEquals(123.36, 0));
  EXPECT_THAT(l9.max(), AlmostEquals(123.54, 0));

  ApproximateQuantity<double> const three = 3.0_⑶;
  EXPECT_THAT(three.min(), AlmostEquals(2.7, 0));
  EXPECT_THAT(three.max(), AlmostEquals(3.3, 0));

  ApproximateQuantity<double> const quote = 11'972_⑴;
  EXPECT_THAT(quote.min(), AlmostEquals(11971, 0));
  EXPECT_THAT(quote.max(), AlmostEquals(11973, 0));
}

TEST(ApproximateQuantityTest, Literals_🄐_🄕) {
  ApproximateQuantity<double> const la = 0x1E3.45p0_🄐;
  EXPECT_THAT(la.min(), AlmostEquals(0x1E3.3Bp0, 0));
  EXPECT_THAT(la.max(), AlmostEquals(0x1E3.4Fp0, 0));

  ApproximateQuantity<double> const lb = 0x1E3.45p0_🄑;
  EXPECT_THAT(lb.min(), AlmostEquals(0x1E3.3Ap0, 0));
  EXPECT_THAT(lb.max(), AlmostEquals(0x1E3.50p0, 0));

  ApproximateQuantity<double> const lc = 0x1E3.45p0_🄒;
  EXPECT_THAT(lc.min(), AlmostEquals(0x1E3.39p0, 0));
  EXPECT_THAT(lc.max(), AlmostEquals(0x1E3.51p0, 0));

  ApproximateQuantity<double> const ld = 0x1E3.45p0_🄓;
  EXPECT_THAT(ld.min(), AlmostEquals(0x1E3.38p0, 0));
  EXPECT_THAT(ld.max(), AlmostEquals(0x1E3.52p0, 0));

  ApproximateQuantity<double> const le = 0x1E3.45p0_🄔;
  EXPECT_THAT(le.min(), AlmostEquals(0x1E3.37p0, 0));
  EXPECT_THAT(le.max(), AlmostEquals(0x1E3.53p0, 0));

  ApproximateQuantity<double> const lf = 0x1E3.45p0_🄕;
  EXPECT_THAT(lf.min(), AlmostEquals(0x1E3.36p0, 0));
  EXPECT_THAT(lf.max(), AlmostEquals(0x1E3.54p0, 0));
}

TEST(ApproximateQuantityTest, Quantities) {
  ApproximateQuantity<Length> const l1 = 123.45_⑴ * Metre;
  EXPECT_THAT(l1.min(), AlmostEquals(123.44 * Metre, 0));
  EXPECT_THAT(l1.max(), AlmostEquals(123.46 * Metre, 1));

  ApproximateQuantity<Area> const l2 = (123.45_⑴ * Metre) * Metre;
  EXPECT_THAT(l2.min(), AlmostEquals(123.44 * Metre * Metre, 0));
  EXPECT_THAT(l2.max(), AlmostEquals(123.46 * Metre * Metre, 1));

  ApproximateQuantity<Frequency> const l3 = 123.45_⑴ / Second;
  EXPECT_THAT(l3.min(), AlmostEquals(123.44 / Second, 0));
  EXPECT_THAT(l3.max(), AlmostEquals(123.46 / Second, 1));

  ApproximateQuantity<Speed> const l4 = 123.45_⑴ * Metre / Second;
  EXPECT_THAT(l4.min(), AlmostEquals(123.44 * Metre / Second, 0));
  EXPECT_THAT(l4.max(), AlmostEquals(123.46 * Metre / Second, 1));
}

TEST(ApproximateQuantityTest, Signs) {
  ApproximateQuantity<double> const l1 = +123.45_⑴;
  EXPECT_THAT(l1.min(), AlmostEquals(123.44, 0));
  EXPECT_THAT(l1.max(), AlmostEquals(123.46, 1));

  ApproximateQuantity<double> const l2 = -123.45_⑴;
  EXPECT_THAT(l2.min(), AlmostEquals(-123.46, 1));
  EXPECT_THAT(l2.max(), AlmostEquals(-123.44, 0));

  ApproximateQuantity<Length> const l3 = +(123.45_⑴ * Metre);
  EXPECT_THAT(l3.min(), AlmostEquals(123.44 * Metre, 0));
  EXPECT_THAT(l3.max(), AlmostEquals(123.46 * Metre, 1));

  ApproximateQuantity<Length> const l4 = -(123.45_⑴ * Metre);
  EXPECT_THAT(l4.min(), AlmostEquals(-123.46 * Metre, 1));
  EXPECT_THAT(l4.max(), AlmostEquals(-123.44 * Metre, 0));
}

TEST(ApproximateQuantityTest, Unit) {
  ApproximateQuantity<Length> const l1 = 123.45_⑴ * Metre;
  EXPECT_EQ(Metre, l1.unit());
  EXPECT_TRUE(l1.has_trivial_unit());

  ApproximateQuantity<Length> const l2 = 123.45_⑴ * (2 * Metre);
  EXPECT_EQ(2 * Metre, l2.unit());
  EXPECT_FALSE(l2.has_trivial_unit());
}

TEST(ApproximateQuantityTest, UlpDistance) {
  ApproximateQuantity<double> const l1 = 123.45_⑴;
  EXPECT_THAT(l1.UlpDistance(123.50), AlmostEquals(5, 3200));

  ApproximateQuantity<double> const l2 = 123.45_⑵;
  EXPECT_THAT(l2.UlpDistance(123.50), AlmostEquals(5, 800));

  ApproximateQuantity<Length> const l3 = 123.45_⑴ * Metre;
  EXPECT_THAT(l3.UlpDistance(123.50 * Metre), AlmostEquals(5, 3200));

  ApproximateQuantity<Length> const l4 = 123.45_⑵ * Metre;
  EXPECT_THAT(l4.UlpDistance(123.50 * Metre), AlmostEquals(5, 800));
}

TEST(ApproximateQuantityTest, DebugString) {
  EXPECT_EQ("123.45(1)", (123.45_⑴).DebugString());
  EXPECT_EQ("123.45(1) m", (123.45_⑴ * Metre).DebugString());
  EXPECT_EQ("123.45(1) * +2.00000000000000000e+00 m",
            (123.45_⑴ * (2 * Metre)).DebugString());

  EXPECT_EQ("-123.45(1)", (-123.45_⑴).DebugString());
  EXPECT_EQ("-123.45(1) m", (-123.45_⑴ * Metre).DebugString());
  EXPECT_EQ("-123.45(1) * +2.00000000000000000e+00 m",
            (-123.45_⑴ * (2 * Metre)).DebugString());
}

}  // namespace testing_utilities
}  // namespace principia
