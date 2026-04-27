#pragma once

/// @file
/// @brief Matrix validation helpers.

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace mathcore
{
    /// Checks whether a matrix is symmetric positive definite within configurable tolerances.
    /// @param matrix Matrix to validate.
    /// @param equality_tolerance Tolerance used for the symmetry comparison.
    /// @param eigenvalue_tolerance Lower bound accepted for the real part of each eigenvalue.
    /// @return `true` when `matrix` is square, approximately symmetric, and all eigenvalues exceed the threshold.
    template <typename Derived>
    [[nodiscard]] inline bool is_symmetric_positive_definite(const Eigen::MatrixBase<Derived> &matrix,
                                                             const typename Derived::Scalar equality_tolerance = typename Derived::Scalar(1e-12),
                                                             const typename Derived::Scalar eigenvalue_tolerance = typename Derived::Scalar(-1e-8))
    {
        if (matrix.rows() != matrix.cols())
        {
            return false;
        }

        const bool is_symmetric = matrix.isApprox(matrix.transpose(), equality_tolerance);
        if (!is_symmetric)
        {
            return false;
        }

        const auto eigenvalues = matrix.eigenvalues();
        return (eigenvalues.real().array() > eigenvalue_tolerance).all();
    }
} // namespace mathcore
