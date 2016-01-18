﻿
#include <sstream>

#include "geometry/r3_element.hpp"
#include "gmock/gmock-matchers.h"
#include "gtest/gtest-death-test.h"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/massive_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/parser.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {

using testing_utilities::AlmostEquals;

namespace quantities {

using si::AstronomicalUnit;
using si::Day;
using si::Degree;
using si::Kilo;
using si::Metre;
using si::Radian;
using si::Second;

class ParserTest : public ::testing::Test {
};

using ParserDeathTest = ParserTest;

TEST_F(ParserTest, SpacesSuccess) {
  EXPECT_EQ(1.23, ParseQuantity<double>("  1.23"));
  EXPECT_EQ(1.23, ParseQuantity<double>("  1.23   "));
  EXPECT_EQ(1.23, ParseQuantity<double>("1.23 "));
  EXPECT_EQ(1.23, ParseQuantity<double>("1.23"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("  1.23 m/s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("  1.23 m/s   "));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m/s "));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m/s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m /s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m/ s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m / s"));
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23m/s"));
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m ^ 3/s^2"));
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m^3 / s^2"));
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m ^ 3 / s ^ 2"));
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m^ 3/s ^2"));
}

TEST_F(ParserDeathTest, SpacesError) {
  EXPECT_DEATH({
    ParseQuantity<double>("1. 23");
  }, "empty");
  EXPECT_DEATH({
    ParseQuantity<double>("1 .23");
  }, "empty");
  EXPECT_DEATH({
    ParseQuantity<Speed>("1. 23 m/s");
  }, "Unsupported.*length");
  EXPECT_DEATH({
    ParseQuantity<Speed>("1 .23 m/s");
  }, "Unsupported.*length");
}

TEST_F(ParserDeathTest, UnitError) {
  EXPECT_DEATH({
    ParseQuantity<Length>("1.23 nm");
  }, "Unsupported.*length");
  EXPECT_DEATH({
    ParseQuantity<Time>("1.23 hr");
  }, "Unsupported.*time");
  EXPECT_DEATH({
    ParseQuantity<Angle>("1.23 grd");
  }, "Unsupported.*angle");
}

TEST_F(ParserTest, ParseDouble) {
  EXPECT_EQ(1.23, ParseQuantity<double>("1.23"));
  EXPECT_EQ(-3.45, ParseQuantity<double>("-3.45"));
}

TEST_F(ParserTest, ParseLength) {
  EXPECT_EQ(1.23 * Metre, ParseQuantity<Length>("1.23 m"));
  EXPECT_EQ(1.23 * Kilo(Metre), ParseQuantity<Length>("1.23 km"));
  EXPECT_EQ(1.23 * AstronomicalUnit, ParseQuantity<Length>("1.23 au"));
}

TEST_F(ParserTest, ParseAngle) {
  EXPECT_EQ(1.23 * Degree, ParseQuantity<Angle>("1.23 deg"));
  EXPECT_EQ(1.23 * Degree, ParseQuantity<Angle>(u8"1.23 °"));
  EXPECT_EQ(1.23 * Radian, ParseQuantity<Angle>("1.23 rad"));
}

TEST_F(ParserTest, ParseSpeed) {
  EXPECT_EQ(1.23 * Metre / Second, ParseQuantity<Speed>("1.23 m/s"));
  EXPECT_EQ(1.23 * Kilo(Metre) / Second, ParseQuantity<Speed>("1.23 km/s"));
  EXPECT_EQ(1.23 * Kilo(Metre) / Day, ParseQuantity<Speed>("1.23 km/d"));
  EXPECT_EQ(1.23 * AstronomicalUnit / Day, ParseQuantity<Speed>("1.23 au/d"));
}

TEST_F(ParserTest, ParseGravitationalParameter) {
  EXPECT_EQ(1.23 * Pow<3>(Metre) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 m^3/s^2"));
  EXPECT_EQ(1.23 * Pow<3>(Kilo(Metre)) / Pow<2>(Second),
            ParseQuantity<GravitationalParameter>("1.23 km^3/s^2"));
  EXPECT_EQ(1.23 * Pow<3>(Kilo(Metre)) / Pow<2>(Day),
            ParseQuantity<GravitationalParameter>("1.23 km^3/d^2"));
  EXPECT_THAT(ParseQuantity<GravitationalParameter>("1.23 au^3/d^2"),
              AlmostEquals(1.23 * Pow<3>(AstronomicalUnit) / Pow<2>(Day), 1));
}

}  // namespace quantities
}  // namespace principia
