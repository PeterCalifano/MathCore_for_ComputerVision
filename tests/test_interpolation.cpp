#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <initializer_list>

#include <mathcore/interpolation/chebyshev.h>

namespace
{
    auto makeYawQuaternion(const double angle) -> Eigen::Vector4d
    {
        const double half_angle = angle / 2.0;

        Eigen::Vector4d quaternion;
        quaternion << std::cos(half_angle), 0.0, 0.0, std::sin(half_angle);
        return quaternion;
    }

    auto makeQuaternionSamples(const Eigen::VectorXd &domain) -> Eigen::MatrixXd
    {
        const double total_angle = 0.8 * std::acos(-1.0);

        Eigen::MatrixXd data(4, domain.size());
        for (Eigen::Index index = 0; index < domain.size(); ++index)
        {
            data.col(index) = makeYawQuaternion(total_angle * domain(index));
        }

        return data;
    }

    void flipQuaternionSigns(Eigen::MatrixXd &data_matrix, const std::initializer_list<Eigen::Index> &columns)
    {
        for (const auto column : columns)
        {
            data_matrix.col(column) *= -1.0;
        }
    }

    struct quaternion_patch_probe : mathcore::CChebyshevQuaternionInterpolator<Eigen::Vector4d>
    {
        using base_type = mathcore::CChebyshevQuaternionInterpolator<Eigen::Vector4d>;
        using matrix_type = base_type::matrix_type;
        using vector_type = base_type::vector_type;

        quaternion_patch_probe(const int polynomial_degree, const vector_type &interpolation_domain)
            : base_type(polynomial_degree, interpolation_domain)
        {
        }

        void applyPatch(matrix_type &data_matrix)
        {
            this->patchDiscontinuities(data_matrix);
        }
    };
} // namespace

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

    mathcore::CChebyshevQuaternionInterpolator<Eigen::Vector4d> interpolator(2, domain, data);
    const auto quaternion = interpolator.evaluateQuaternion(0.5);

    REQUIRE(quaternion.norm() == Catch::Approx(1.0).margin(1e-12));
    REQUIRE(quaternion(3) == Catch::Approx(std::sqrt(0.5)).margin(1e-8));
}

TEST_CASE("patchDiscontinuities_removes_consecutive_quaternion_sign_flips", "[interpolation]")
{
    const Eigen::VectorXd domain = Eigen::VectorXd::LinSpaced(5, 0.0, 1.0);
    const Eigen::MatrixXd expected = makeQuaternionSamples(domain);

    Eigen::MatrixXd sign_switched = expected;
    flipQuaternionSigns(sign_switched, {1, 2, 4});

    quaternion_patch_probe probe(4, domain);
    probe.applyPatch(sign_switched);

    REQUIRE(sign_switched.isApprox(expected, 1e-12));

    for (Eigen::Index index = 1; index < sign_switched.cols(); ++index)
    {
        REQUIRE(sign_switched.col(index - 1).dot(sign_switched.col(index)) > 0.0);
    }
}

TEST_CASE("CChebyshevQuaternionInterpolator_applies_sign_switch_patching_when_enabled", "[interpolation]")
{
    const Eigen::VectorXd domain = Eigen::VectorXd::LinSpaced(6, 0.0, 1.0);
    const Eigen::MatrixXd reference_data = makeQuaternionSamples(domain);

    Eigen::MatrixXd sign_switched_data = reference_data;
    flipQuaternionSigns(sign_switched_data, {1, 2, 4});

    mathcore::CChebyshevQuaternionInterpolator<Eigen::Vector4d> reference(5, domain, reference_data);
    mathcore::CChebyshevQuaternionInterpolator<Eigen::Vector4d> patched(5, domain, true, sign_switched_data);

    for (const double evaluation_time : std::array<double, 3>{0.15, 0.50, 0.85})
    {
        const auto expected = reference.evaluateQuaternion(evaluation_time);
        const auto actual = patched.evaluateQuaternion(evaluation_time);

        REQUIRE(mathcore::quaternionDistance(actual, expected) < 1e-12);
    }
}

TEST_CASE("CChebyshevQuaternionInterpolator_refit_applies_sign_switch_patching_by_default", "[interpolation]")
{
    const Eigen::VectorXd domain = Eigen::VectorXd::LinSpaced(6, 0.0, 1.0);
    const Eigen::MatrixXd reference_data = makeQuaternionSamples(domain);

    Eigen::MatrixXd sign_switched_data = reference_data;
    flipQuaternionSigns(sign_switched_data, {1, 2, 4});

    mathcore::CChebyshevQuaternionInterpolator<Eigen::Vector4d> reference(5, domain, reference_data);
    mathcore::CChebyshevQuaternionInterpolator<Eigen::Vector4d> refit_interpolator(5, domain);
    refit_interpolator.fitCoefficientMatrix(sign_switched_data, domain);

    for (const double evaluation_time : std::array<double, 3>{0.15, 0.50, 0.85})
    {
        const auto expected = reference.evaluateQuaternion(evaluation_time);
        const auto actual = refit_interpolator.evaluateQuaternion(evaluation_time);

        REQUIRE(mathcore::quaternionDistance(actual, expected) < 1e-12);
    }
}

TEST_CASE("CChebyshevQuaternionInterpolator_throws_for_empty_refit_domain", "[interpolation]")
{
    Eigen::VectorXd domain(2);
    domain << 0.0, 1.0;

    Eigen::MatrixXd data(4, 2);
    data.col(0) << 1.0, 0.0, 0.0, 0.0;
    data.col(1) << 0.0, 0.0, 0.0, 1.0;

    mathcore::CChebyshevQuaternionInterpolator<Eigen::Vector4d> interpolator(1, domain, data);

    Eigen::MatrixXd empty_data(4, 0);
    const Eigen::VectorXd empty_domain;

    REQUIRE_THROWS_AS(interpolator.fitCoefficientMatrix(empty_data, empty_domain), std::invalid_argument);
}
