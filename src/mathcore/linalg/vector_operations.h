#pragma once

/// @file
/// @brief Header defining templates of common vector operations.

#include <mathcore/linalg/vector_traits.h>

namespace mathcore
{
    /// Builds the skew-symmetric matrix associated with a 3D vector.
    /// @param vector Input 3D vector.
    /// @return The matrix `S(v)` such that `S(v) * w == v x w`.
    template <vector3_like Vector>
    [[nodiscard]] inline auto skew_symmetric_matrix(const Vector &vector)
        -> Eigen::Matrix<vector_scalar_t<Vector>, 3, 3>
    {
        // TODO optimize by avoiding the temporary Eigen vector and directly accessing the coefficients through the traits
        const auto eigen_vector = to_eigen_vector3(vector);
        using Scalar = vector_scalar_t<Vector>;

        Eigen::Matrix<Scalar, 3, 3> result;
        result << Scalar(0), -eigen_vector(2), eigen_vector(1),
            eigen_vector(2), Scalar(0), -eigen_vector(0),
            -eigen_vector(1), eigen_vector(0), Scalar(0);

        return result;
    }
} // namespace mathcore
