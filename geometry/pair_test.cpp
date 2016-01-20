﻿
#include <sstream>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/identity.hpp"
#include "geometry/pair.hpp"
#include "geometry/point.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "gtest/gtest-death-test.h"
#include "gtest/gtest.h"
#include "integrators/ordinary_differential_equations.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "serialization/geometry.pb.h"

// Define this to check that the illegal operations are actually rejected.
#undef CHECK_ILLEGAL

namespace principia {

using quantities::Action;
using quantities::Angle;
using quantities::Dimensions;
using quantities::Energy;
using quantities::Entropy;
using quantities::Frequency;
using quantities::Quantity;
using quantities::SolidAngle;
using quantities::Time;
using quantities::Winding;

namespace geometry {

class PairTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;

  struct Universe;

  PairTest()
      : p1_(P1() + V1({4 * SIUnit<Action>(),
                       5 * SIUnit<Action>(),
                       6 * SIUnit<Action>()})),
        p2_(P2() + V2({16 * SIUnit<Winding>(),
                       15 * SIUnit<Winding>(),
                       14 * SIUnit<Winding>()})),
        v1_({1 * SIUnit<Action>(),
             2 * SIUnit<Action>(),
             3 * SIUnit<Action>()}),
        v2_({13 * SIUnit<Winding>(),
             12 * SIUnit<Winding>(),
             11 * SIUnit<Winding>()}),
        pp_(p1_, p2_),
        pv_(p1_, v2_),
        vp_(v1_, p2_),
        vv_(v1_, v2_) {}

  using V1 = Vector<Action, World>;
  using P1 = Point<V1>;
  using V2 = Vector<Winding, World>;
  using P2 = Point<V2>;
  using PP = Pair<P1, P2>;
  using PV = Pair<P1, V2>;
  using VP = Pair<V1, P2>;
  using VV = Pair<V1, V2>;

  P1 p1_;
  P2 p2_;
  V1 v1_;
  V2 v2_;
  PP pp_;
  PV pv_;
  VP vp_;
  VV vv_;
};

using PairDeathTest = PairTest;

TEST_F(PairTest, MemberAddition) {
  EXPECT_EQ(PP(P1() + V1({5 * SIUnit<Action>(),
                          7 * SIUnit<Action>(),
                          9 * SIUnit<Action>()}),
               P2() + V2({29 * SIUnit<Winding>(),
                          27 * SIUnit<Winding>(),
                          25 * SIUnit<Winding>()})),
            pp_ + vv_);
  EXPECT_EQ(PV(P1() + V1({5 * SIUnit<Action>(),
                          7 * SIUnit<Action>(),
                          9 * SIUnit<Action>()}),
               V2({26 * SIUnit<Winding>(),
                   24 * SIUnit<Winding>(),
                   22 * SIUnit<Winding>()})),
            pv_ + vv_);
  EXPECT_EQ(VP(V1({2 * SIUnit<Action>(),
                   4 * SIUnit<Action>(),
                   6 * SIUnit<Action>()}),
               P2() + V2({29 * SIUnit<Winding>(),
                          27 * SIUnit<Winding>(),
                          25 * SIUnit<Winding>()})),
            vp_ + vv_);
  EXPECT_EQ(VV(V1({2 * SIUnit<Action>(),
                   4 * SIUnit<Action>(),
                   6 * SIUnit<Action>()}),
               V2({26 * SIUnit<Winding>(),
                   24 * SIUnit<Winding>(),
                   22 * SIUnit<Winding>()})),
            vv_ + vv_);
}

TEST_F(PairTest, MemberSubtraction) {
  EXPECT_EQ(PP(P1() + V1({3 * SIUnit<Action>(),
                          3 * SIUnit<Action>(),
                          3 * SIUnit<Action>()}),
               P2() + V2({3 * SIUnit<Winding>(),
                          3 * SIUnit<Winding>(),
                          3 * SIUnit<Winding>()})),
            pp_ - vv_);
  EXPECT_EQ(PV(P1() + V1({3 * SIUnit<Action>(),
                          3 * SIUnit<Action>(),
                          3 * SIUnit<Action>()}),
               V2({0 * SIUnit<Winding>(),
                   0 * SIUnit<Winding>(),
                   0 * SIUnit<Winding>()})),
            pv_ - vv_);
  EXPECT_EQ(VP(V1({0 * SIUnit<Action>(),
                   0 * SIUnit<Action>(),
                   0 * SIUnit<Action>()}),
               P2() + V2({3 * SIUnit<Winding>(),
                          3 * SIUnit<Winding>(),
                          3 * SIUnit<Winding>()})),
            vp_ - vv_);
  EXPECT_EQ(VV(V1({0 * SIUnit<Action>(),
                   0 * SIUnit<Action>(),
                   0 * SIUnit<Action>()}),
               V2({0 * SIUnit<Winding>(),
                   0 * SIUnit<Winding>(),
                   0 * SIUnit<Winding>()})),
            vv_ - vv_);
}

TEST_F(PairTest, MemberAdditionTo) {
  pp_ += vv_;
  EXPECT_EQ(PP(P1() + V1({5 * SIUnit<Action>(),
                          7 * SIUnit<Action>(),
                          9 * SIUnit<Action>()}),
               P2() + V2({29 * SIUnit<Winding>(),
                          27 * SIUnit<Winding>(),
                          25 * SIUnit<Winding>()})),
            pp_);
  pv_ += vv_;
  EXPECT_EQ(PV(P1() + V1({5 * SIUnit<Action>(),
                          7 * SIUnit<Action>(),
                          9 * SIUnit<Action>()}),
               V2({26 * SIUnit<Winding>(),
                   24 * SIUnit<Winding>(),
                   22 * SIUnit<Winding>()})),
            pv_);
  vp_ += vv_;
  EXPECT_EQ(VP(V1({2 * SIUnit<Action>(),
                   4 * SIUnit<Action>(),
                   6 * SIUnit<Action>()}),
               P2() + V2({29 * SIUnit<Winding>(),
                          27 * SIUnit<Winding>(),
                          25 * SIUnit<Winding>()})),
            vp_);
  vv_ += vv_;
  EXPECT_EQ(VV(V1({2 * SIUnit<Action>(),
                   4 * SIUnit<Action>(),
                   6 * SIUnit<Action>()}),
               V2({26 * SIUnit<Winding>(),
                   24 * SIUnit<Winding>(),
                   22 * SIUnit<Winding>()})),
            vv_);
}

TEST_F(PairTest, MemberSubtractionFrom) {
  pp_ -= vv_;
  EXPECT_EQ(PP(P1() + V1({3 * SIUnit<Action>(),
                          3 * SIUnit<Action>(),
                          3 * SIUnit<Action>()}),
               P2() + V2({3 * SIUnit<Winding>(),
                          3 * SIUnit<Winding>(),
                          3 * SIUnit<Winding>()})),
            pp_);
  pv_ -= vv_;
  EXPECT_EQ(PV(P1() + V1({3 * SIUnit<Action>(),
                          3 * SIUnit<Action>(),
                          3 * SIUnit<Action>()}),
               V2({0 * SIUnit<Winding>(),
                   0 * SIUnit<Winding>(),
                   0 * SIUnit<Winding>()})),
            pv_);
  vp_ -= vv_;
  EXPECT_EQ(VP(V1({0 * SIUnit<Action>(),
                   0 * SIUnit<Action>(),
                   0 * SIUnit<Action>()}),
               P2() + V2({3 * SIUnit<Winding>(),
                          3 * SIUnit<Winding>(),
                          3 * SIUnit<Winding>()})),
            vp_);
  vv_ -= vv_;
  EXPECT_EQ(VV(V1({0 * SIUnit<Action>(),
                   0 * SIUnit<Action>(),
                   0 * SIUnit<Action>()}),
               V2({0 * SIUnit<Winding>(),
                   0 * SIUnit<Winding>(),
                   0 * SIUnit<Winding>()})),
            vv_);
}

TEST_F(PairTest, MemberEquality) {
  EXPECT_TRUE(pp_ == PP(P1() + V1({4 * SIUnit<Action>(),
                                   5 * SIUnit<Action>(),
                                   6 * SIUnit<Action>()}),
                        P2() + V2({16 * SIUnit<Winding>(),
                                   15 * SIUnit<Winding>(),
                                   14 * SIUnit<Winding>()})));
  EXPECT_FALSE(pp_ == PP(P1() + V1({4 * SIUnit<Action>(),
                                    5 * SIUnit<Action>(),
                                    6 * SIUnit<Action>()}),
                         P2() + V2({26 * SIUnit<Winding>(),
                                15 * SIUnit<Winding>(),
                                14 * SIUnit<Winding>()})));
  EXPECT_TRUE(pv_ == PV(P1() + V1({4 * SIUnit<Action>(),
                                   5 * SIUnit<Action>(),
                                   6 * SIUnit<Action>()}),
                        V2({13 * SIUnit<Winding>(),
                            12 * SIUnit<Winding>(),
                            11 * SIUnit<Winding>()})));
  EXPECT_FALSE(pv_ == PV(P1() + V1({4 * SIUnit<Action>(),
                                    15 * SIUnit<Action>(),
                                    6 * SIUnit<Action>()}),
                         V2({13 * SIUnit<Winding>(),
                             12 * SIUnit<Winding>(),
                             11 * SIUnit<Winding>()})));
  EXPECT_TRUE(vp_ == VP(V1({1 * SIUnit<Action>(),
                            2 * SIUnit<Action>(),
                            3 * SIUnit<Action>()}),
                        P2() + V2({16 * SIUnit<Winding>(),
                               15 * SIUnit<Winding>(),
                               14 * SIUnit<Winding>()})));
  EXPECT_FALSE(vp_ == VP(V1({1 * SIUnit<Action>(),
                             2 * SIUnit<Action>(),
                             13 * SIUnit<Action>()}),
                         P2() + V2({16 * SIUnit<Winding>(),
                                15 * SIUnit<Winding>(),
                                14 * SIUnit<Winding>()})));
  EXPECT_TRUE(vv_ == VV(V1({1 * SIUnit<Action>(),
                            2 * SIUnit<Action>(),
                            3 * SIUnit<Action>()}),
                        V2({13 * SIUnit<Winding>(),
                            12 * SIUnit<Winding>(),
                            11 * SIUnit<Winding>()})));
  EXPECT_FALSE(vv_ == VV(V1({1 * SIUnit<Action>(),
                             2 * SIUnit<Action>(),
                             3 * SIUnit<Action>()}),
                         V2({13 * SIUnit<Winding>(),
                             12 * SIUnit<Winding>(),
                             1 * SIUnit<Winding>()})));
}

TEST_F(PairTest, MemberInequality) {
  EXPECT_FALSE(pp_ != PP(P1() + V1({4 * SIUnit<Action>(),
                                    5 * SIUnit<Action>(),
                                    6 * SIUnit<Action>()}),
                         P2() + V2({16 * SIUnit<Winding>(),
                                    15 * SIUnit<Winding>(),
                                    14 * SIUnit<Winding>()})));
  EXPECT_TRUE(pp_ != PP(P1() + V1({4 * SIUnit<Action>(),
                                   5 * SIUnit<Action>(),
                                   6 * SIUnit<Action>()}),
                        P2() + V2({26 * SIUnit<Winding>(),
                                   15 * SIUnit<Winding>(),
                                   14 * SIUnit<Winding>()})));
  EXPECT_FALSE(pv_ != PV(P1() + V1({4 * SIUnit<Action>(),
                                    5 * SIUnit<Action>(),
                                    6 * SIUnit<Action>()}),
                         V2({13 * SIUnit<Winding>(),
                             12 * SIUnit<Winding>(),
                             11 * SIUnit<Winding>()})));
  EXPECT_TRUE(pv_ != PV(P1() + V1({4 * SIUnit<Action>(),
                                   15 * SIUnit<Action>(),
                                   6 * SIUnit<Action>()}),
                        V2({13 * SIUnit<Winding>(),
                            12 * SIUnit<Winding>(),
                            11 * SIUnit<Winding>()})));
  EXPECT_FALSE(vp_ != VP(V1({1 * SIUnit<Action>(),
                             2 * SIUnit<Action>(),
                             3 * SIUnit<Action>()}),
                         P2() + V2({16 * SIUnit<Winding>(),
                                15 * SIUnit<Winding>(),
                                14 * SIUnit<Winding>()})));
  EXPECT_TRUE(vp_ != VP(V1({1 * SIUnit<Action>(),
                            2 * SIUnit<Action>(),
                            13 * SIUnit<Action>()}),
                        P2() + V2({16 * SIUnit<Winding>(),
                               15 * SIUnit<Winding>(),
                               14 * SIUnit<Winding>()})));
  EXPECT_FALSE(vv_ != VV(V1({1 * SIUnit<Action>(),
                             2 * SIUnit<Action>(),
                             3 * SIUnit<Action>()}),
                         V2({13 * SIUnit<Winding>(),
                             12 * SIUnit<Winding>(),
                             11 * SIUnit<Winding>()})));
  EXPECT_TRUE(vv_ != VV(V1({1 * SIUnit<Action>(),
                            2 * SIUnit<Action>(),
                            3 * SIUnit<Action>()}),
                        V2({13 * SIUnit<Winding>(),
                            12 * SIUnit<Winding>(),
                            1 * SIUnit<Winding>()})));
}

TEST_F(PairTest, AffineSubtraction) {
  PP const pp = pp_ + vv_;
  EXPECT_EQ(vv_, pp - pp_);
  PV const pv = pv_ + vv_;
  EXPECT_EQ(vv_, pv - pv_);
  VP const vp = vp_ + vv_;
  EXPECT_EQ(vv_, vp - vp_);
  // No test for VV, that would be a vector subtraction.
}

TEST_F(PairTest, UnaryPlus) {
  VV const vv = +vv_;
  EXPECT_EQ(vv_, vv);
#ifdef CHECK_ILLEGAL
  auto const pp = +pp_;
  auto const pv = +pv_;
  auto const vp = +vp_;
#endif
}

TEST_F(PairTest, UnaryMinus) {
  VV const vv = -vv_;
  EXPECT_EQ(VV(V1({-1 * SIUnit<Action>(),
                   -2 * SIUnit<Action>(),
                   -3 * SIUnit<Action>()}),
               V2({-13 * SIUnit<Winding>(),
                   -12 * SIUnit<Winding>(),
                   -11 * SIUnit<Winding>()})),
            vv);
#ifdef CHECK_ILLEGAL
  auto const pp = -pp_;
  auto const pv = -pv_;
  auto const vp = -vp_;
#endif
}

TEST_F(PairTest, LeftMultiplication) {
  EXPECT_EQ(VV(V1({3 * SIUnit<Action>(),
                   6 * SIUnit<Action>(),
                   9 * SIUnit<Action>()}),
               V2({39 * SIUnit<Winding>(),
                   36 * SIUnit<Winding>(),
                   33 * SIUnit<Winding>()})),
            3 * vv_);
#ifdef CHECK_ILLEGAL
  auto const pp = 3 * pp_;
  auto const pv = 3 * pv_;
  auto const vp = 3 * vp_;
#endif
}

TEST_F(PairTest, RightMultiplication) {
  using Whatever1 = Quantity<Dimensions<2, 1, -1, 0, 0, 0, 0, 0, 0, 1>>;
  using Whatever2 = Quantity<Dimensions<0, 0, 0, 0, 0, 0, 0, 1, 0, 1>>;
  using VWhatever1 = Vector<Whatever1, World>;
  using VWhatever2 = Vector<Whatever2, World>;
  using Pear = Pair<VWhatever1, VWhatever2>;
  EXPECT_EQ(Pear(VWhatever1({1.5 * SIUnit<Whatever1>(),
                             3.0 * SIUnit<Whatever1>(),
                             4.5 * SIUnit<Whatever1>()}),
                 VWhatever2({19.5 * SIUnit<Whatever2>(),
                             18.0 * SIUnit<Whatever2>(),
                             16.5 * SIUnit<Whatever2>()})),
            vv_ * (1.5 * SIUnit<SolidAngle>()));
#ifdef CHECK_ILLEGAL
  auto const pp = pp_ * (1.5 * SIUnit<SolidAngle>());
  auto const pv = pv_ * (1.5 * SIUnit<SolidAngle>());
  auto const vp = vp_ * (1.5 * SIUnit<SolidAngle>());
#endif
}

TEST_F(PairTest, RightDivision) {
  using VEnergy = Vector<Energy, World>;
  using VFrequency = Vector<Frequency, World>;
  using Pear = Pair<VEnergy, VFrequency>;
  EXPECT_EQ(Pear(VEnergy({0.5 * SIUnit<Energy>(),
                          1.0 * SIUnit<Energy>(),
                          1.5 * SIUnit<Energy>()}),
                 VFrequency({6.5 * SIUnit<Frequency>(),
                             6.0 * SIUnit<Frequency>(),
                             5.5 * SIUnit<Frequency>()})),
            vv_ / (2 * SIUnit<Time>()));
#ifdef CHECK_ILLEGAL
  auto const pp = pp_ / (2 * SIUnit<Time>());
  auto const pv = pv_ / (2 * SIUnit<Time>());
  auto const vp = vp_ / (2 * SIUnit<Time>());
#endif
}

TEST_F(PairTest, MultiplicationBy) {
  vv_ *= 4;
  EXPECT_EQ(VV(V1({4 * SIUnit<Action>(),
                   8 * SIUnit<Action>(),
                   12 * SIUnit<Action>()}),
               V2({52 * SIUnit<Winding>(),
                   48 * SIUnit<Winding>(),
                   44 * SIUnit<Winding>()})),
            vv_);
#ifdef CHECK_ILLEGAL
  pp_ *= 4;
  pv_ *= 4;
  vp_ *= 4;
#endif
}

TEST_F(PairTest, DivisionBy) {
  vv_ /= 0.25;
  EXPECT_EQ(VV(V1({4 * SIUnit<Action>(),
                   8 * SIUnit<Action>(),
                   12 * SIUnit<Action>()}),
               V2({52 * SIUnit<Winding>(),
                   48 * SIUnit<Winding>(),
                   44 * SIUnit<Winding>()})),
            vv_);
#ifdef CHECK_ILLEGAL
  pp_ /= 0.25;
  pv_ /= 0.25;
  vp_ /= 0.25;
#endif
}

TEST_F(PairTest, Streaming) {
  LOG(ERROR) << "pp_ = " << pp_;
  LOG(ERROR) << "pv_ = " << pv_;
  LOG(ERROR) << "vp_ = " << vp_;
  LOG(ERROR) << "vv_ = " << vv_;
}

TEST_F(PairDeathTest, SerializationError) {
  serialization::Pair message;
  pv_.WriteToMessage(&message);
  EXPECT_DEATH({
    PP const pp = PP::ReadFromMessage(message);
  }, "has_point");
  EXPECT_DEATH({
    VP const vp = VP::ReadFromMessage(message);
  }, "has_multivector");
  EXPECT_DEATH({
    VV const vv = VV::ReadFromMessage(message);
  }, "has_multivector");
}

TEST_F(PairTest, SerializationSuccess) {
  serialization::Pair message;

  pp_.WriteToMessage(&message);
  EXPECT_TRUE(message.t1().has_point());
  EXPECT_FALSE(message.t1().has_multivector());
  EXPECT_TRUE(message.t2().has_point());
  EXPECT_FALSE(message.t2().has_multivector());
  PP const pp = PP::ReadFromMessage(message);
  EXPECT_EQ(pp_, pp);

  pv_.WriteToMessage(&message);
  EXPECT_TRUE(message.t1().has_point());
  EXPECT_FALSE(message.t1().has_multivector());
  EXPECT_FALSE(message.t2().has_point());
  EXPECT_TRUE(message.t2().has_multivector());
  PV const pv = PV::ReadFromMessage(message);
  EXPECT_EQ(pv_, pv);

  vp_.WriteToMessage(&message);
  EXPECT_FALSE(message.t1().has_point());
  EXPECT_TRUE(message.t1().has_multivector());
  EXPECT_TRUE(message.t2().has_point());
  EXPECT_FALSE(message.t2().has_multivector());
  VP const vp = VP::ReadFromMessage(message);
  EXPECT_EQ(vp_, vp);

  vv_.WriteToMessage(&message);
  EXPECT_FALSE(message.t1().has_point());
  EXPECT_TRUE(message.t1().has_multivector());
  EXPECT_FALSE(message.t2().has_point());
  EXPECT_TRUE(message.t2().has_multivector());
  VV const vv = VV::ReadFromMessage(message);
  EXPECT_EQ(vv_, vv);
}

TEST_F(PairDeathTest, BarycentreCalculatorError) {
  using PPBarycentreCalculator = BarycentreCalculator<PP, Entropy>;
  using PVBarycentreCalculator = BarycentreCalculator<PV, Entropy>;
  using VPBarycentreCalculator = BarycentreCalculator<VP, Entropy>;
  using VVBarycentreCalculator = BarycentreCalculator<VV, Entropy>;
  EXPECT_DEATH({
    PPBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
  EXPECT_DEATH({
    PVBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
  EXPECT_DEATH({
    VPBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
  EXPECT_DEATH({
    VVBarycentreCalculator calculator;
    calculator.Get();
  }, "Empty BarycentreCalculator");
}

// Because we can't access individual members this test doesn't fully exercise
// the computations, so we'll redo some testing for DegreesOfFreedom.
TEST_F(PairTest, BarycentreCalculatorSuccess) {
  {
    BarycentreCalculator<PP, double> calculator;
    calculator.Add(pp_, 3);
    PP barycentre = calculator.Get();
    EXPECT_EQ(pp_, barycentre);
    calculator.Add(pp_ + vv_, 5);
    barycentre = calculator.Get();
    EXPECT_EQ(pp_ + 5.0 / 8.0 * vv_, barycentre);
    calculator.Add(pp_ - 3 * vv_, 8);
    barycentre = calculator.Get();
    EXPECT_EQ(pp_ - 19.0 / 16.0 * vv_, barycentre);
  }
  {
    BarycentreCalculator<PV, double> calculator;
    calculator.Add(pv_, 3);
    PV barycentre = calculator.Get();
    EXPECT_EQ(pv_, barycentre);
    calculator.Add(pv_ + vv_, 5);
    barycentre = calculator.Get();
    EXPECT_EQ(pv_ + 5.0 / 8.0 * vv_, barycentre);
    calculator.Add(pv_ - 3 * vv_, 8);
    barycentre = calculator.Get();
    EXPECT_EQ(pv_ - 19.0 / 16.0 * vv_, barycentre);
  }
  {
    BarycentreCalculator<VP, double> calculator;
    calculator.Add(vp_, 3);
    VP barycentre = calculator.Get();
    EXPECT_EQ(vp_, barycentre);
    calculator.Add(vp_ + vv_, 5);
    barycentre = calculator.Get();
    EXPECT_EQ(vp_ + 5.0 / 8.0 * vv_, barycentre);
    calculator.Add(vp_ - 3 * vv_, 8);
    barycentre = calculator.Get();
    EXPECT_EQ(vp_ - 19.0 / 16.0 * vv_, barycentre);
  }
  {
    BarycentreCalculator<VV, double> calculator;
    calculator.Add(vv_, 3);
    VV barycentre = calculator.Get();
    EXPECT_EQ(vv_, barycentre);
    calculator.Add(vv_ + vv_, 5);
    barycentre = calculator.Get();
    EXPECT_EQ(13.0 / 8.0 * vv_, barycentre);
    calculator.Add(-2 * vv_, 8);
    barycentre = calculator.Get();
    EXPECT_EQ(-3.0 / 16.0 * vv_, barycentre);
  }
}

TEST_F(PairTest, Mappable) {
  using VV1U = Vector<Action, Universe>;
  using VV2U = Vector<Winding, Universe>;
  using VVU = Pair<VV1U, VV2U>;
  VVU const vv = Identity<World, Universe>()(vv_);
  EXPECT_EQ(VVU(VV1U({1 * SIUnit<Action>(),
                      2 * SIUnit<Action>(),
                      3 * SIUnit<Action>()}),
                VV2U({13 * SIUnit<Winding>(),
                      12 * SIUnit<Winding>(),
                      11 * SIUnit<Winding>()})),
            vv);

#ifdef CHECK_ILLEGAL
  auto const pp = Identity<World, Universe>()(pp_);
  auto const pv = Identity<World, Universe>()(pv_);
  auto const vp = Identity<World, Universe>()(vp_);
#endif
}

}  // namespace geometry
}  // namespace principia
