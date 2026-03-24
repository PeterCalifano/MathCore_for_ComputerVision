#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <mathcore/random/scatter.h>

TEST_CASE("scatter_vector3_returns_mean_when_stddev_is_non_positive", "[random]")
{
    const Eigen::Vector3d mean{1.0, 2.0, 3.0};
    const auto sample = mathcore::scatter_vector3<Eigen::Vector3d>(0.0, mean, 7);

    REQUIRE(sample.isApprox(mean, 1e-12));
}

TEST_CASE("scatter_vector3_is_seed_deterministic", "[random]")
{
    const Eigen::Vector3d mean{0.0, 0.0, 0.0};
    const auto lhs = mathcore::scatter_vector3<Eigen::Vector3d>(0.1, mean, 42);
    const auto rhs = mathcore::scatter_vector3<Eigen::Vector3d>(0.1, mean, 42);

    REQUIRE(lhs.isApprox(rhs, 1e-12));
}
