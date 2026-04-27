#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <mathcore/linalg/skew_symmetric.h>
#include <mathcore/random/vector_sampling.h>
#include <mathcore/rotations/quaternion.h>

namespace
{
struct custom_vec3
{
    double x{};
    double y{};
    double z{};
};

struct custom_quaternion
{
    double x{};
    double y{};
    double z{};
    double w{};
};
} // namespace

template <>
struct mathcore::vector_traits<custom_vec3>
{
    using scalar_type = double;
    static constexpr bool is_specialized = true;
    static constexpr std::size_t static_size = 3;

    static auto get(const custom_vec3 &value, const std::size_t index) -> scalar_type
    {
        switch (index)
        {
        case 0:
            return value.x;
        case 1:
            return value.y;
        default:
            return value.z;
        }
    }

    static auto make(const Eigen::Matrix<scalar_type, 3, 1> &coefficients) -> custom_vec3
    {
        return {coefficients(0), coefficients(1), coefficients(2)};
    }

    static auto zero() -> custom_vec3
    {
        return {};
    }
};

template <>
struct mathcore::SQuatTraits<custom_quaternion>
{
    using scalar_type = double;
    using convention = mathcore::SScalarLast;

    static auto toScalarFirst(const custom_quaternion &value) -> Eigen::Matrix<scalar_type, 4, 1>
    {
        Eigen::Matrix<scalar_type, 4, 1> coefficients;
        coefficients << value.w, value.x, value.y, value.z;
        return coefficients;
    }

    static auto fromScalarFirst(const Eigen::Matrix<scalar_type, 4, 1> &coefficients) -> custom_quaternion
    {
        return {coefficients(1), coefficients(2), coefficients(3), coefficients(0)};
    }
};

TEST_CASE("custom_vector_type_can_use_linalg_and_random_helpers", "[extensibility]")
{
    const custom_vec3 vector{1.0, 2.0, 3.0};
    const auto skew = mathcore::skew_symmetric_matrix(vector);
    const auto scattered = mathcore::scatter_vector3<custom_vec3>(0.0, vector, 7);

    REQUIRE(skew(1, 0) == Catch::Approx(3.0));
    REQUIRE(scattered.x == Catch::Approx(1.0));
    REQUIRE(scattered.y == Catch::Approx(2.0));
}

TEST_CASE("custom_quaternion_type_can_use_rotation_helpers", "[extensibility]")
{
    const custom_quaternion identity{0.0, 0.0, 0.0, 1.0};
    const auto halfway = mathcore::slerp(identity, identity, 0.5);
    const auto dcm = mathcore::quaternionToDcm(identity);

    REQUIRE(halfway.w == Catch::Approx(1.0));
    REQUIRE(dcm.isApprox(Eigen::Matrix3d::Identity(), 1e-12));
}
