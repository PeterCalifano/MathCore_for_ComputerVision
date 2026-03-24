/**
 * @file interpolator.h
 * @author PeterC (petercalifano.gs@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-03-24
 */
#pragma once

/// @file
/// @brief Base utilities for polynomial interpolation over sampled data.

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace mathcore
{
    /// Randomly samples matched points from a data matrix and interpolation domain.
    /// @param num_test_points Maximum number of points to sample.
    /// @param data_points Output sampled columns from `data_matrix`.
    /// @param test_points Output sampled times from `interpolation_domain`.
    /// @param data_matrix Matrix whose columns are the samples to be selected.
    /// @param interpolation_domain Domain values associated with `data_matrix`.
    template <std::floating_point Scalar>
    inline void sampleTestPoints(const int numTestPoints,
                                 std::vector<Eigen::VectorX<Scalar>> &dataPoints,
                                 Eigen::VectorX<Scalar> &testPoints,
                                 const Eigen::MatrixX<Scalar> &dataMatrix,
                                 const Eigen::VectorX<Scalar> &interpolationDomain)
    {
        if (dataMatrix.cols() != interpolationDomain.size())
        {
            throw std::invalid_argument("Data matrix and interpolation domain must have matching lengths.");
        }

        if (dataMatrix.cols() == 0)
        {
            dataPoints.clear();
            testPoints.resize(0);
            return;
        }

        const int effectiveNumTestPoints = std::min(numTestPoints, static_cast<int>(interpolationDomain.size()));

        // Generate a random permutation of indices corresponding to the interpolation domain samples
        std::vector<int> allIndices(static_cast<std::size_t>(interpolationDomain.size()));
        std::iota(allIndices.begin(), allIndices.end(), 0);

        std::random_device device;
        std::mt19937 generator(device());
        std::shuffle(allIndices.begin(), allIndices.end(), generator);

        // Allocate samples arrays
        dataPoints.resize(static_cast<std::size_t>(effectiveNumTestPoints));
        testPoints.resize(effectiveNumTestPoints);

        // Extract the randomly selected samples
        for (int pointIndex = 0; pointIndex < effectiveNumTestPoints; ++pointIndex)
        {
            const auto sampledIndex = allIndices[static_cast<std::size_t>(pointIndex)];
            dataPoints.at(static_cast<std::size_t>(pointIndex)) = dataMatrix.col(sampledIndex); // ALlocate entry in std::vector of Eigen vectors
            testPoints(pointIndex) = interpolationDomain(sampledIndex);                         // Allocate entry in Eigen vector
        }
    }

    // CLASS IMPLEMENTATIONS
    /**
     * @brief Base class for polynomial interpolators over sampled vector data.
     *
     * @tparam Scalar
     */
    template <std::floating_point Scalar = double>
    class CInterpolatorBase
    {
        // TYPEDEFS
      public:
        using scalar_type = Scalar;
        using vector_type = Eigen::VectorX<Scalar>;
        using matrix_type = Eigen::MatrixX<Scalar>;

        // CONSTRUCTORS
        CInterpolatorBase() = default;

        /// Constructs an interpolator for a fixed domain and output dimensionality.
        CInterpolatorBase(const int polynomialDegree,
                          const vector_type &interpolationDomain,
                          const int num_components)
            : interpolation_domain_(interpolationDomain),
              polynomial_degree_(polynomialDegree),
              num_components_(num_components)
        {
            domain_size_ = static_cast<int>(interpolation_domain_.size());
            if (domain_size_ < polynomial_degree_ + 1)
            {
                throw std::invalid_argument("Interpolation domain must be at least polynomial_degree + 1 samples long.");
            }

            time_lower_bound_ = interpolation_domain_.minCoeff();
            time_upper_bound_ = interpolation_domain_.maxCoeff();
            scaled_interpolation_domain_ = scaleInterpolationDomain(interpolation_domain_);
        }

        // DESTRUCTOR
        virtual ~CInterpolatorBase() = default;

        // SETTERS

        /// Enables or disables automatic post-fit validation
        void setAutoCheck(const bool auto_check)
        {
            auto_check_ = auto_check;
        }

        // GETTERS

        /// Returns the fitted coefficient matrix.
        [[nodiscard]] auto coefficients() const -> matrix_type
        {
            return coefficient_matrix_;
        }

        /// Returns the lower bound of the interpolation domain.
        [[nodiscard]] auto timeLowerBound() const -> scalar_type
        {
            return time_lower_bound_;
        }

        /// Returns the upper bound of the interpolation domain.
        [[nodiscard]] auto timeUpperBound() const -> scalar_type
        {
            return time_upper_bound_;
        }

        /// Returns the configured polynomial degree.
        [[nodiscard]] auto polynomialDegree() const -> int
        {
            return polynomial_degree_;
        }

        /// Returns the number of output components per sample.
        [[nodiscard]] auto numComponents() const -> int
        {
            return num_components_;
        }

        /// Returns the number of points in the interpolation domain.
        [[nodiscard]] auto domainSize() const -> int
        {
            return domain_size_;
        }

        /// Returns the original interpolation domain.
        [[nodiscard]] auto interpolationDomain() const -> vector_type
        {
            return interpolation_domain_;
        }

        /// Returns the interpolation domain scaled to `[-1, 1]`.
        [[nodiscard]] auto scaledInterpolationDomain() const -> vector_type
        {
            return scaled_interpolation_domain_;
        }

        /// Indicates whether interpolation coefficients have been fitted.
        [[nodiscard]] auto hasCoefficients() const -> bool
        {
            return coefficient_matrix_.size() != 0;
        }

        // METHODS
        /**
         * @brief Fits the coefficient matrix for the supplied samples and domain.
         *
         * @param data_matrix
         * @param interpolation_domain
         */
        virtual void fitCoefficientMatrix(const matrix_type &data_matrix, const vector_type &interpolation_domain) = 0;

        /**
         * @brief Evaluates the interpolant at a point in the domain
         *
         * @param eval_point
         * @param scale_value
         * @return vector_type
         */
        [[nodiscard]] virtual vector_type evaluate(const scalar_type eval_point, const bool scale_value = true) const = 0;

        /**
         * @brief Evaluates the polynomial basis at a point in the scaled domain
         *
         * @param eval_point
         * @return vector_type
         */
        [[nodiscard]] virtual vector_type evaluatePolynomialBase(const scalar_type eval_point) const = 0;

        /// Compares interpolation output against reference samples.
        /// @param data_points Ground-truth vectors to compare against.
        /// @param evaluation_points Domain locations associated with `data_points`.
        /// @param absolute_tolerance Maximum allowed absolute error norm.
        /// @param relative_tolerance Maximum allowed relative error norm.
        /// @param scale_value Whether `evaluation_points` should be mapped to the normalized interpolation interval.
        /// @return `true` if either the absolute or relative error criterion passes.
        [[nodiscard]] bool checkInterpolationError(const std::vector<vector_type> &data_points,
                                                   const vector_type &evaluation_points,
                                                   const scalar_type absolute_tolerance = scalar_type(1e-5),
                                                   const scalar_type relative_tolerance = scalar_type(1e-4),
                                                   const bool scale_value = true) const
        {
            if (!hasCoefficients())
            {
                throw std::runtime_error("Polynomial coefficient matrix is empty.");
            }

            if (static_cast<int>(data_points.size()) != evaluation_points.size())
            {
                throw std::invalid_argument("Data points and evaluation points must have the same length.");
            }

            vector_type absolute_errors = vector_type::Zero(evaluation_points.size());
            vector_type relative_errors = vector_type::Zero(evaluation_points.size());

            // Compute errors for each evaluation point
            for (Eigen::Index index = 0; index < evaluation_points.size(); ++index)
            {
                // Absolute error
                const auto interpolated = evaluate(evaluation_points(index), scale_value);
                absolute_errors(index) = (data_points.at(static_cast<std::size_t>(index)) - interpolated).norm();

                // Relative error
                const auto reference_norm = data_points.at(static_cast<std::size_t>(index)).norm();
                if (reference_norm > scalar_type(10) * std::numeric_limits<scalar_type>::epsilon())
                {
                    relative_errors(index) = absolute_errors(index) / reference_norm;
                }
            }

            // Determine validity
            const bool absolute_condition = absolute_errors.maxCoeff() < absolute_tolerance;
            const bool relative_condition = relative_errors.maxCoeff() <= scalar_type(0) ||
                                            relative_errors.maxCoeff() < relative_tolerance;
            return absolute_condition || relative_condition;
        }

        /**
         * @brief Removes sign flips from a quaternion sample sequence before fitting
         *
         * @param quaternion_matrix
         */
        void patchDiscontinuities(Eigen::Matrix<Scalar, 4, Eigen::Dynamic> &quaternion_matrix)
        {
            if (quaternion_matrix.rows() != 4)
            {
                throw std::invalid_argument("Quaternion matrix must be 4xN.");
            }

            const std::size_t num_samples = static_cast<std::size_t>(quaternion_matrix.cols());
            if (num_samples == 0)
            {
                switch_intervals_.clear();
                num_switches_ = 0;
                return;
            }

            std::vector<bool> sign_switch_detection_mask(num_samples, false);
            Eigen::Matrix<Scalar, 4, 1> previous = quaternion_matrix.col(0);

            // TODO: review implementation and test more

            // Loop through quaternion samples and detect sign flips by checking the dot product between consecutive quaternions
            for (std::size_t column = 1; column < num_samples; ++column)
            {
                Eigen::Matrix<Scalar, 4, 1> current = quaternion_matrix.col(static_cast<Eigen::Index>(column));

                // If the dot product between the current and previous quaternion is negative, sign flip
                if (previous.dot(current) < scalar_type(0))
                {
                    current = -current;
                    sign_switch_detection_mask[column] = true;
                }

                quaternion_matrix.col(static_cast<Eigen::Index>(column)) = current;
                previous = current;
            }

            num_switches_ = static_cast<std::uint8_t>(
                std::count(sign_switch_detection_mask.begin(), sign_switch_detection_mask.end(), true));

            switch_intervals_.clear();
            if (num_switches_ == 0)
            {
                return;
            }

            // Identify contiguous intervals of sign flips
            std::vector<std::size_t> transitions;
            for (std::size_t column = 1; column < num_samples; ++column)
            {
                if (sign_switch_detection_mask[column] != sign_switch_detection_mask[column - 1])
                {
                    transitions.push_back(column);
                }
            }

            if (sign_switch_detection_mask.front())
            {
                transitions.insert(transitions.begin(), 0);
            }

            if (sign_switch_detection_mask.back())
            {
                transitions.push_back(num_samples);
            }

            for (std::size_t index = 0; index + 1 < transitions.size(); index += 2)
            {
                switch_intervals_.emplace_back(transitions[index], transitions[index + 1] - 1);
            }
        }

        /// Scales a domain from its native range into the Chebyshev interval `[-1, 1]`.
        [[nodiscard]] static auto scaleInterpolationDomain(const vector_type &interpolation_domain) -> vector_type
        {
            if (interpolation_domain.size() == 0)
            {
                return {};
            }

            const scalar_type lower_bound = interpolation_domain.minCoeff();
            const scalar_type upper_bound = interpolation_domain.maxCoeff();
            if (std::abs(upper_bound - lower_bound) <= std::numeric_limits<scalar_type>::epsilon())
            {
                return vector_type::Zero(interpolation_domain.size());
            }

            vector_type scaled_domain(interpolation_domain.size());
            for (Eigen::Index index = 0; index < interpolation_domain.size(); ++index)
            {
                scaled_domain(index) =
                    (scalar_type(2) * interpolation_domain(index) - (upper_bound + lower_bound)) /
                    (upper_bound - lower_bound);
            }

            return scaled_domain;
        }

        /// Tests whether a scaled evaluation point falls inside a detected switch interval.
        /// @param scaled_eval_point Evaluation point expressed in the scaled domain.
        /// @param scaled_interpolation_domain Sample times scaled to `[-1, 1]`.
        /// @param switch_intervals Closed index intervals where sign changes were detected.
        /// @return `true` if the point maps into one of the intervals.
        [[nodiscard]] static bool isEvalPointInSwitchInterval(const scalar_type scaled_eval_point,
                                                              const vector_type &scaled_interpolation_domain,
                                                              const std::vector<std::pair<std::size_t, std::size_t>> &switch_intervals)
        {
            if (scaled_interpolation_domain.size() == 0)
            {
                return false;
            }

            if (scaled_eval_point < scaled_interpolation_domain(0) ||
                scaled_eval_point > scaled_interpolation_domain(scaled_interpolation_domain.size() - 1))
            {
                return false;
            }

            const auto *data = scaled_interpolation_domain.data();
            const auto num_times = static_cast<std::size_t>(scaled_interpolation_domain.size());
            const auto *iterator = std::lower_bound(data, data + num_times, scaled_eval_point);
            const std::size_t index = std::min<std::size_t>(static_cast<std::size_t>(std::distance(data, iterator)),
                                                            num_times - 1);

            return std::ranges::any_of(switch_intervals,
                                       [index](const auto &interval)
                                       {
                                           return index >= interval.first && index <= interval.second;
                                       });
        }

      protected:
        // PROTECTED METHODS
        /// Checks whether an evaluation point lies inside the configured domain bounds.
        [[nodiscard]] bool checkEvalPointValidity(const scalar_type eval_point) const
        {
            return eval_point >= time_lower_bound_ && eval_point <= time_upper_bound_;
        }

        /// Validates that the data matrix shape matches the interpolator configuration.
        [[nodiscard]] bool checkDataMatrixValidity(const matrix_type &data_matrix) const
        {
            if (num_components_ != data_matrix.rows())
            {
                throw std::invalid_argument("Data matrix row count does not match the configured number of components.");
            }

            if (data_matrix.cols() != interpolation_domain_.size())
            {
                throw std::invalid_argument("Data matrix and interpolation domain lengths do not match.");
            }

            if (data_matrix.cols() < polynomial_degree_ + 1)
            {
                throw std::invalid_argument("Data matrix must contain at least polynomial_degree + 1 samples.");
            }

            return true;
        }

      protected:
        // PROTECTED DATA MEMBERS
        int domain_size_{0};
        vector_type interpolation_domain_{};
        vector_type scaled_interpolation_domain_{};
        scalar_type time_lower_bound_{scalar_type(0)};
        scalar_type time_upper_bound_{scalar_type(0)};
        int polynomial_degree_{0};
        int num_components_{0};
        bool auto_check_{true};
        matrix_type coefficient_matrix_{};
        std::vector<std::pair<std::size_t, std::size_t>> switch_intervals_{};
        std::uint8_t num_switches_{0};
    };
} // namespace mathcore
