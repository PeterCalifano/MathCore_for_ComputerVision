#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <array>
#include <cmath>

#include <mathcore/rotations/quaternion.h>

TEST_CASE("dcm_round_trip_preserves_identity_rotation", "[rotations]")
{
    const auto quaternion = mathcore::dcmToQuaternion<Eigen::Vector4d>(Eigen::Matrix3d::Identity());
    const auto dcm = mathcore::quaternionToDcm(quaternion);

    REQUIRE(quaternion(0) == Catch::Approx(1.0));
    REQUIRE(quaternion(1) == Catch::Approx(0.0));
    REQUIRE(dcm.isApprox(Eigen::Matrix3d::Identity(), 1e-12));
}

TEST_CASE("slerp_halfway_matches_expected_quaternion", "[rotations]")
{
    Eigen::Vector4d start;
    start << 1.0, 0.0, 0.0, 0.0;

    const double angle = std::acos(-1.0);
    Eigen::Vector4d target;
    target << std::cos(angle / 2.0), 0.0, 0.0, std::sin(angle / 2.0);

    const auto interpolated = mathcore::slerp(start, target, 0.5);

    REQUIRE(interpolated(0) == Catch::Approx(std::sqrt(0.5)).margin(1e-12));
    REQUIRE(interpolated(3) == Catch::Approx(std::sqrt(0.5)).margin(1e-12));
}

TEST_CASE("quaternion_distance_is_zero_for_equivalent_sign_flip", "[rotations]")
{
    Eigen::Vector4d lhs;
    lhs << 1.0, 0.0, 0.0, 0.0;
    const Eigen::Vector4d rhs = -lhs;

    REQUIRE(mathcore::quaternionDistance(lhs, rhs) == Catch::Approx(0.0));
}

TEST_CASE("std_array_quaternion_support_round_trips_and_interpolates", "[rotations]")
{
    const std::array<double, 4> identity{1.0, 0.0, 0.0, 0.0};
    const auto dcm = mathcore::quaternionToDcm(identity);
    const auto round_trip = mathcore::dcmToQuaternion<std::array<double, 4>>(dcm);
    const auto halfway = mathcore::slerp(identity, identity, 0.5);

    REQUIRE(dcm.isApprox(Eigen::Matrix3d::Identity(), 1e-12));
    REQUIRE(round_trip[0] == Catch::Approx(1.0));
    REQUIRE(round_trip[1] == Catch::Approx(0.0));
    REQUIRE(halfway[0] == Catch::Approx(1.0));
    REQUIRE(mathcore::quaternionDistance(identity, round_trip) == Catch::Approx(0.0));
}
