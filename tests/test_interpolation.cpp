#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>

#include <mathcore/interpolation/chebyshev.h>

TEST_CASE("CChebyshevInterpolator_fits_linear_signal", "[interpolation]")
{
    Eigen::VectorXd domain(5);
    domain << -1.0, -0.5, 0.0, 0.5, 1.0;

    Eigen::MatrixXd data(1, 5);
    data << -1.0, -0.5, 0.0, 0.5, 1.0;

    mathcore::CChebyshevInterpolator<double> interpolator(1, domain, 1, false);
    interpolator.fitCoefficientMatrix(data, domain);

    const auto value = interpolator.evaluate(0.25);
    REQUIRE(value(0) == Catch::Approx(0.25).margin(1e-10));
}

TEST_CASE("CChebyshevQuaternionInterpolator_returns_normalized_quaternion", "[interpolation]")
{
    Eigen::VectorXd domain(3);
    domain << 0.0, 0.5, 1.0;

    Eigen::MatrixXd data(4, 3);
    data.col(0) << 1.0, 0.0, 0.0, 0.0;
    data.col(1) << std::sqrt(0.5), 0.0, 0.0, std::sqrt(0.5);
    data.col(2) << 0.0, 0.0, 0.0, 1.0;

    mathcore::CChebyshevQuaternionInterpolator<Eigen::Vector4d> interpolator(2, domain, true, data);
    const auto quaternion = interpolator.evaluateQuaternion(0.5);

    REQUIRE(quaternion.norm() == Catch::Approx(1.0).margin(1e-12));
    REQUIRE(quaternion(3) == Catch::Approx(std::sqrt(0.5)).margin(1e-8));
}
