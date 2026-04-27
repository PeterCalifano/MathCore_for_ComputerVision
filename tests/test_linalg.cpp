#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <mathcore/linalg/vector_operations.h>
#include <mathcore/linalg/validation.h>

TEST_CASE("skew_symmetric_matrix_matches_cross_product_form", "[linalg]")
{
    const Eigen::Vector3d vector{1.0, 2.0, 3.0};
    const auto skew = mathcore::skew_symmetric_matrix(vector);

    REQUIRE(skew(0, 1) == Catch::Approx(-3.0));
    REQUIRE(skew(0, 2) == Catch::Approx(2.0));
    REQUIRE(skew(2, 1) == Catch::Approx(1.0));
}

TEST_CASE("spd_validation_accepts_spd_matrix", "[linalg]")
{
    const Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity();
    REQUIRE(mathcore::is_symmetric_positive_definite(matrix));
}
