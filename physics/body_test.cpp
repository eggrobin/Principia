#include "physics/body.hpp"

#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/oblate_body.hpp"
#include "physics/rotating_body.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"

namespace principia {

using geometry::AngularVelocity;
using geometry::Frame;
using geometry::Normalize;
using si::Radian;
using si::Second;
using ::testing::IsNull;
using ::testing::NotNull;

namespace physics {

class BodyTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST, true>;

  // We need that so the comma doesn't get caught in macros.
  using Direction = Vector<double, World>;

  template<typename Tag, Tag tag>
  void TestOblateBody() {
    using F = Frame<Tag, tag, true>;

    AngularVelocity<F> const angular_velocity =
        AngularVelocity<F>({-1 * Radian / Second,
                            2 * Radian / Second,
                            5 * Radian / Second});
    auto const oblate_body =
        OblateBody<F>(17 * SIUnit<GravitationalParameter>(),
                      RotatingBody<F>::Parameters(1 * Radian,
                                                  Instant(),
                                                  angular_velocity),
                      OblateBody<F>::Parameters(
                          163 * SIUnit<Order2ZonalCoefficient>()));

    serialization::Body message;
    OblateBody<F> const* cast_oblate_body;
    oblate_body.WriteToMessage(&message);
    EXPECT_TRUE(message.has_massive_body());
    EXPECT_FALSE(message.has_massless_body());
    EXPECT_TRUE(message.massive_body().HasExtension(
                    serialization::OblateBody::oblate_body));

    not_null<std::unique_ptr<MassiveBody const>> const massive_body =
        MassiveBody::ReadFromMessage(message);
    EXPECT_EQ(oblate_body.gravitational_parameter(),
              massive_body->gravitational_parameter());
    cast_oblate_body = dynamic_cast<OblateBody<F> const*>(&*massive_body);
    EXPECT_THAT(cast_oblate_body, NotNull());
  }

  AngularVelocity<World> angular_velocity_ =
      AngularVelocity<World>({-1 * Radian / Second,
                              2 * Radian / Second,
                              5 * Radian / Second});
  MasslessBody massless_body_;
  MassiveBody massive_body_ =
      MassiveBody(42 * SIUnit<GravitationalParameter>());
  OblateBody<World> oblate_body_ =
      OblateBody<World>(17 * SIUnit<GravitationalParameter>(),
                        RotatingBody<World>::Parameters(1 * Radian,
                                                        Instant(),
                                                        angular_velocity_),
                        OblateBody<World>::Parameters(
                            163 * SIUnit<Order2ZonalCoefficient>()));
};

using BodyDeathTest = BodyTest;

TEST_F(BodyTest, MasslessSerializationSuccess) {
  EXPECT_TRUE(massless_body_.is_massless());
  EXPECT_FALSE(massless_body_.is_oblate());

  serialization::Body message;
  MasslessBody const* cast_massless_body;
  massless_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massless_body());
  EXPECT_FALSE(message.has_massive_body());

  // Direct deserialization.
  // No members to test in this class, we just check that it doesn't crash.
  massless_body_ = *MasslessBody::ReadFromMessage(message);

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  // NOTE(egg): The &* is a quick way to explicitly forget |not_null|ness. We
  // cannot strip the |not_null| from the previous line because MSVC does not
  // support move conversion at the moment.
  cast_massless_body = dynamic_cast<MasslessBody const*>(&*body);
  EXPECT_THAT(cast_massless_body, NotNull());
}

// The best serialization revenge.
TEST_F(BodyTest, MassiveSerializationSuccess) {
  EXPECT_FALSE(massive_body_.is_massless());
  EXPECT_FALSE(massive_body_.is_oblate());

  serialization::Body message;
  MassiveBody const* cast_massive_body;
  massive_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massive_body());
  EXPECT_FALSE(message.has_massless_body());
  EXPECT_EQ(42, message.massive_body().gravitational_parameter().magnitude());

  // Direct deserialization.
  MassiveBody const massive_body = *MassiveBody::ReadFromMessage(message);
  EXPECT_EQ(massive_body_.gravitational_parameter(),
            massive_body.gravitational_parameter());

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body>> body = Body::ReadFromMessage(message);
  cast_massive_body = dynamic_cast<MassiveBody*>(&*body);
  EXPECT_THAT(cast_massive_body, NotNull());
  EXPECT_EQ(massive_body_.gravitational_parameter(),
            cast_massive_body->gravitational_parameter());
}

TEST_F(BodyTest, OblateSerializationSuccess) {
  EXPECT_FALSE(oblate_body_.is_massless());
  EXPECT_TRUE(oblate_body_.is_oblate());

  serialization::Body message;
  OblateBody<World> const* cast_oblate_body;
  oblate_body_.WriteToMessage(&message);
  EXPECT_TRUE(message.has_massive_body());
  EXPECT_FALSE(message.has_massless_body());
  EXPECT_TRUE(
      message.massive_body().HasExtension(
          serialization::OblateBody::oblate_body));
  EXPECT_EQ(17, message.massive_body().gravitational_parameter().magnitude());
  serialization::OblateBody const oblateness_information =
      message.massive_body().GetExtension(
          serialization::OblateBody::oblate_body);
  EXPECT_EQ(163, oblateness_information.j2().magnitude());
  EXPECT_EQ(Normalize(angular_velocity_),
            Direction::ReadFromMessage(oblateness_information.axis()));

  // Direct deserialization.
  OblateBody<World> const oblate_body =
      *OblateBody<World>::ReadFromMessage(message);
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            oblate_body.gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), oblate_body.j2());
  EXPECT_EQ(oblate_body_.axis(), oblate_body.axis());

  // Dispatching from |MassiveBody|.
  not_null<std::unique_ptr<MassiveBody const>> const massive_body =
      MassiveBody::ReadFromMessage(message);
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            massive_body->gravitational_parameter());
  cast_oblate_body = dynamic_cast<OblateBody<World> const*>(&*massive_body);
  EXPECT_THAT(cast_oblate_body, NotNull());
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            cast_oblate_body->gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), cast_oblate_body->j2());
  EXPECT_EQ(oblate_body_.axis(), cast_oblate_body->axis());

  // Dispatching from |Body|.
  not_null<std::unique_ptr<Body const>> const body =
      Body::ReadFromMessage(message);
  cast_oblate_body = dynamic_cast<OblateBody<World> const*>(&*body);
  EXPECT_THAT(cast_oblate_body, NotNull());
  EXPECT_EQ(oblate_body_.gravitational_parameter(),
            cast_oblate_body->gravitational_parameter());
  EXPECT_EQ(oblate_body_.j2(), cast_oblate_body->j2());
  EXPECT_EQ(oblate_body_.axis(), cast_oblate_body->axis());
}

TEST_F(BodyTest, AllFrames) {
  TestOblateBody<serialization::Frame::PluginTag,
                 serialization::Frame::ALICE_SUN>();
  TestOblateBody<serialization::Frame::PluginTag,
                 serialization::Frame::ALICE_WORLD>();
  TestOblateBody<serialization::Frame::PluginTag,
                 serialization::Frame::BARYCENTRIC>();
  TestOblateBody<serialization::Frame::PluginTag,
                 serialization::Frame::PRE_BOREL_BARYCENTRIC>();
  TestOblateBody<serialization::Frame::PluginTag,
                 serialization::Frame::RENDERING>();
  TestOblateBody<serialization::Frame::PluginTag,
                 serialization::Frame::WORLD>();
  TestOblateBody<serialization::Frame::PluginTag,
                 serialization::Frame::WORLD_SUN>();

  TestOblateBody<serialization::Frame::SolarSystemTag,
                 serialization::Frame::ICRF_J2000_ECLIPTIC>();
  TestOblateBody<serialization::Frame::SolarSystemTag,
                 serialization::Frame::ICRF_J2000_EQUATOR>();

  TestOblateBody<serialization::Frame::TestTag,
                 serialization::Frame::TEST>();
  TestOblateBody<serialization::Frame::TestTag,
                 serialization::Frame::TEST1>();
  TestOblateBody<serialization::Frame::TestTag,
                 serialization::Frame::TEST2>();
  TestOblateBody<serialization::Frame::TestTag,
                 serialization::Frame::FROM>();
  TestOblateBody<serialization::Frame::TestTag,
                 serialization::Frame::THROUGH>();
  TestOblateBody<serialization::Frame::TestTag,
                 serialization::Frame::TO>();
}

}  // namespace physics
}  // namespace principia
