#pragma once

/// @file
/// @brief Quaternion distance, conversion, composition, and interpolation helpers.

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include <mathcore/rotations/quaternion_traits.h>

namespace mathcore
{
    /// Computes a sign-invariant distance measure between two quaternions.
    /// @param lhs First quaternion.
    /// @param rhs Second quaternion.
    /// @return `1 - |dot(lhs, rhs)|`, which is zero for equal orientations.
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto quaternionDistance(const Quaternion &lhs, const Quaternion &rhs) -> quaternion_scalar_t<Quaternion>
    {
        const auto lhs_coefficients = getCoefficientsScalarFirst(lhs);
        const auto rhs_coefficients = getCoefficientsScalarFirst(rhs);
        return quaternion_scalar_t<Quaternion>(1) - std::abs(lhs_coefficients.dot(rhs_coefficients));
    }

    /// Alias for `quaternionDistance`.
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto quatDotProductDistance(const Quaternion &lhs, const Quaternion &rhs) -> quaternion_scalar_t<Quaternion>
    {
        return quaternionDistance(lhs, rhs);
    }

    /// Multiplies two quaternions using Hamilton product rules.
    /// @param lhs Left-hand quaternion.
    /// @param rhs Right-hand quaternion.
    /// @return The composed rotation `lhs * rhs`.
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto quaternionProduct(const Quaternion &lhs, const Quaternion &rhs)
        -> std::remove_cvref_t<Quaternion>
    {
        const auto lhs_coefficients = getCoefficientsScalarFirst(lhs);
        const auto rhs_coefficients = getCoefficientsScalarFirst(rhs);
        using Scalar = quaternion_scalar_t<Quaternion>;

        Eigen::Matrix<Scalar, 4, 1> result;
        result(0) = lhs_coefficients(0) * rhs_coefficients(0) - lhs_coefficients(1) * rhs_coefficients(1) -
                    lhs_coefficients(2) * rhs_coefficients(2) - lhs_coefficients(3) * rhs_coefficients(3);
        result(1) = lhs_coefficients(0) * rhs_coefficients(1) + lhs_coefficients(1) * rhs_coefficients(0) +
                    lhs_coefficients(2) * rhs_coefficients(3) - lhs_coefficients(3) * rhs_coefficients(2);
        result(2) = lhs_coefficients(0) * rhs_coefficients(2) - lhs_coefficients(1) * rhs_coefficients(3) +
                    lhs_coefficients(2) * rhs_coefficients(0) + lhs_coefficients(3) * rhs_coefficients(1);
        result(3) = lhs_coefficients(0) * rhs_coefficients(3) + lhs_coefficients(1) * rhs_coefficients(2) -
                    lhs_coefficients(2) * rhs_coefficients(1) + lhs_coefficients(3) * rhs_coefficients(0);

        return makeQuatScalarFirst<Quaternion>(result);
    }

    /// Alias for `quaternionProduct`.
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto quatCross(const Quaternion &lhs, const Quaternion &rhs)
        -> std::remove_cvref_t<Quaternion>
    {
        return quaternionProduct(lhs, rhs);
    }

    /// Converts a quaternion into a 3x3 direction cosine matrix.
    /// @param quaternion Quaternion-like input value.
    /// @return The rotation matrix corresponding to `quaternion`.
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto quaternionToDcm(const Quaternion &quaternion)
        -> Eigen::Matrix<quaternion_scalar_t<Quaternion>, 3, 3>
    {
        const auto coefficients = getCoefficientsScalarFirst(quaternion);
        const auto qs = coefficients(0);
        const auto qx = coefficients(1);
        const auto qy = coefficients(2);
        const auto qz = coefficients(3);

        Eigen::Matrix<quaternion_scalar_t<Quaternion>, 3, 3> dcm;
        dcm(0, 0) = qs * qs + qx * qx - qy * qy - qz * qz;
        dcm(1, 0) = 2 * (qx * qy - qz * qs);
        dcm(2, 0) = 2 * (qx * qz + qy * qs);

        dcm(0, 1) = 2 * (qx * qy + qz * qs);
        dcm(1, 1) = qs * qs - qx * qx + qy * qy - qz * qz;
        dcm(2, 1) = 2 * (qy * qz - qx * qs);

        dcm(0, 2) = 2 * (qx * qz - qy * qs);
        dcm(1, 2) = 2 * (qy * qz + qx * qs);
        dcm(2, 2) = qs * qs - qx * qx - qy * qy + qz * qz;

        return dcm;
    }

    /// Alias for `quaternionToDcm`.
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto quatToDcm(const Quaternion &quaternion)
        -> Eigen::Matrix<quaternion_scalar_t<Quaternion>, 3, 3>
    {
        return quaternionToDcm(quaternion);
    }

    /// Converts a 3x3 direction cosine matrix into a quaternion.
    /// @param matrix Rotation matrix to convert.
    /// @return A quaternion of type `OutputQuaternion` with scalar-first semantics.
    template <TQuatLike OutputQuaternion = Eigen::Vector4d, typename Derived>
    [[nodiscard]] inline auto dcmToQuaternion(const Eigen::MatrixBase<Derived> &matrix)
        -> std::remove_cvref_t<OutputQuaternion>
    {
        static_assert(Derived::RowsAtCompileTime == 3 || Derived::RowsAtCompileTime == Eigen::Dynamic,
                      "dcmToQuaternion expects a 3x3 matrix.");
        static_assert(Derived::ColsAtCompileTime == 3 || Derived::ColsAtCompileTime == Eigen::Dynamic,
                      "dcmToQuaternion expects a 3x3 matrix.");

        using Scalar = quaternion_scalar_t<OutputQuaternion>;
        const Eigen::Matrix<Scalar, 3, 3> dcm = matrix.template cast<Scalar>();
        const Scalar dcm_trace = dcm(0, 0) + dcm(1, 1) + dcm(2, 2);

        std::array<Scalar, 4> candidates{};
        // Select the reconstruction branch from the dominant diagonal/trace term. This avoids
        // dividing by a nearly zero quantity when the rotation is close to a singular case.
        candidates[0] = dcm_trace;
        candidates[1] = dcm(0, 0);
        candidates[2] = dcm(1, 1);
        candidates[3] = dcm(2, 2);

        const int index_of_maximum = static_cast<int>(std::distance(candidates.begin(),
                                                                    std::max_element(candidates.begin(), candidates.end())));

        Eigen::Matrix<Scalar, 4, 1> scalar_last_output;
        // The closed-form branches below reconstruct the component associated with the largest
        // candidate first, then recover the remaining coefficients from the off-diagonal entries.
        switch (index_of_maximum)
        {
        case 1:
            scalar_last_output << Scalar(1) + Scalar(2) * candidates[1] - dcm_trace,
                dcm(0, 1) + dcm(1, 0),
                dcm(0, 2) + dcm(2, 0),
                dcm(1, 2) - dcm(2, 1);
            break;
        case 2:
            scalar_last_output << dcm(1, 0) + dcm(0, 1),
                Scalar(1) + Scalar(2) * candidates[2] - dcm_trace,
                dcm(1, 2) + dcm(2, 1),
                dcm(2, 0) - dcm(0, 2);
            break;
        case 3:
            scalar_last_output << dcm(2, 0) + dcm(0, 2),
                dcm(2, 1) + dcm(1, 2),
                Scalar(1) + Scalar(2) * candidates[3] - dcm_trace,
                dcm(0, 1) - dcm(1, 0);
            break;
        case 0:
            scalar_last_output << dcm(1, 2) - dcm(2, 1),
                dcm(2, 0) - dcm(0, 2),
                dcm(0, 1) - dcm(1, 0),
                Scalar(1) + dcm_trace;
            break;
        default:
            throw std::runtime_error("dcmToQuaternion: invalid DCM input.");
        }

        scalar_last_output.normalize();

        Eigen::Matrix<Scalar, 4, 1> scalar_first_output;
        scalar_first_output << scalar_last_output(3), scalar_last_output(0), scalar_last_output(1), scalar_last_output(2);

        return makeQuatScalarFirst<OutputQuaternion>(scalar_first_output);
    }

    /// Performs spherical linear interpolation between two quaternions.
    /// @param start Quaternion at `normalized_time == 0`.
    /// @param target Quaternion at `normalized_time == 1`.
    /// @param normalized_time Interpolation factor, typically in `[0, 1]`.
    /// @param linear_threshold Dot-product threshold above which normalized linear interpolation is used.
    /// @return The interpolated quaternion.
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto slerp(const Quaternion &start,
                                    const Quaternion &target,
                                    const quaternion_scalar_t<Quaternion> normalized_time,
                                    const quaternion_scalar_t<Quaternion> linear_threshold = quaternion_scalar_t<Quaternion>(0.9995))
        -> std::remove_cvref_t<Quaternion>
    {
        using Scalar = quaternion_scalar_t<Quaternion>;

        auto start_coefficients = getCoefficientsScalarFirst(start);
        auto target_coefficients = getCoefficientsScalarFirst(target);

        Scalar cosine = start_coefficients.dot(target_coefficients);
        if (cosine < Scalar(0))
        {
            // q and -q represent the same orientation; flipping keeps interpolation on the
            // shortest arc of the unit 4D sphere.
            target_coefficients = -target_coefficients;
            cosine = -cosine;
        }

        if (cosine > linear_threshold)
        {
            // Near-coincident endpoints make the spherical basis poorly conditioned, so fall
            // back to normalized lerp and keep the result on the unit sphere.
            const auto interpolated =
                (start_coefficients + normalized_time * (target_coefficients - start_coefficients)).normalized();
            return makeQuatScalarFirst<Quaternion>(interpolated);
        }

        cosine = std::clamp(cosine, Scalar(-1), Scalar(1));

        const Scalar theta_0 = std::acos(cosine);
        const Scalar theta = theta_0 * normalized_time;
        // Remove the component of the target along the start quaternion to obtain the unit
        // tangent direction inside the interpolation plane.
        const auto orthogonal = (target_coefficients - start_coefficients * cosine).normalized();
        const auto interpolated = start_coefficients * std::cos(theta) + orthogonal * std::sin(theta);

        return makeQuatScalarFirst<Quaternion>(interpolated);
    }

    /// Alias for `slerp`.
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto sphericalLerp(const Quaternion &start,
                                            const Quaternion &target,
                                            const quaternion_scalar_t<Quaternion> normalized_time,
                                            const quaternion_scalar_t<Quaternion> linear_threshold = quaternion_scalar_t<Quaternion>(0.9995))
        -> std::remove_cvref_t<Quaternion>
    {
        return slerp(start, target, normalized_time, linear_threshold);
    }
} // namespace mathcore
