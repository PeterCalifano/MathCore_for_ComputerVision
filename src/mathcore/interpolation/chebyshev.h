#pragma once

/// @file
/// @brief Chebyshev-polynomial interpolators for vector and quaternion data.

#include <mathcore/interpolation/interpolator.h>
#include <mathcore/rotations/quaternion.h>

namespace mathcore
{
    /// Fits and evaluates a Chebyshev polynomial interpolant for vector-valued data.
    template <std::floating_point Scalar = double>
    class CChebyshevInterpolator : public CInterpolatorBase<Scalar>
    {
        // TYPEDEFS
      public:
        using base_type = CInterpolatorBase<Scalar>;
        using typename base_type::matrix_type;
        using typename base_type::scalar_type;
        using typename base_type::vector_type;

        // CONSTRUCTORS
        CChebyshevInterpolator() = default;

        /// Constructs an interpolator for the supplied domain and output dimensionality.
        CChebyshevInterpolator(const int polynomial_degree,
                               const vector_type &interpolation_domain,
                               const int num_components,
                               const bool auto_check = true)
            : base_type(polynomial_degree, interpolation_domain, num_components)
        {
            this->setAutoCheck(auto_check);
        }

        /// Constructs an interpolator and immediately fits it to the supplied samples.
        CChebyshevInterpolator(const int polynomial_degree,
                               const vector_type &interpolation_domain,
                               const int num_components,
                               const matrix_type &data_matrix,
                               const bool auto_check = true)
            : base_type(polynomial_degree, interpolation_domain, num_components)
        {
            this->setAutoCheck(auto_check);
            fitCoefficientMatrix(data_matrix, interpolation_domain);
        }

        // METHODS

        /// Fits the Chebyshev coefficient matrix for the currently configured domain.
        void fitCoefficientMatrix(const matrix_type &data_matrix)
        {
            (void)this->checkDataMatrixValidity(data_matrix);

            const auto scaled_domain = base_type::scaleInterpolationDomain(this->interpolation_domain_);
            this->coefficient_matrix_.resize(this->num_components_, this->polynomial_degree_ + 1);

            matrix_type regression_matrix(this->polynomial_degree_ + 1, data_matrix.cols());
            for (Eigen::Index time_index = 0; time_index < data_matrix.cols(); ++time_index)
            {
                // Each column stores the Chebyshev basis evaluated at one sample time. The final
                // coefficient matrix is solved row-wise against this shared basis matrix.
                regression_matrix.col(time_index) = evaluatePolynomialBase(scaled_domain(time_index));
            }

            // Solve regression_matrix^T * coeffs^T = data_matrix^T in least-squares form. SVD is
            // used here because it is robust to ill-conditioned sample layouts.
            Eigen::JacobiSVD<matrix_type> svd(regression_matrix.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
            this->coefficient_matrix_ = svd.solve(data_matrix.transpose()).transpose();

            if (!this->auto_check_)
            {
                return;
            }

            std::vector<vector_type> sampled_data_points;
            vector_type sampled_evaluation_points;
            sampleTestPoints<Scalar>(std::min(10, static_cast<int>(data_matrix.cols())),
                                     sampled_data_points,
                                     sampled_evaluation_points,
                                     data_matrix,
                                     this->interpolation_domain_);

            if (!this->checkInterpolationError(sampled_data_points,
                                               sampled_evaluation_points,
                                               scalar_type(1e-5),
                                               scalar_type(1e-6),
                                               true))
            {
                throw std::runtime_error("Chebyshev interpolation polynomial fitting failed.");
            }
        }

        /// Fits the Chebyshev coefficient matrix after updating the interpolation domain.
        void fitCoefficientMatrix(const matrix_type &data_matrix, const vector_type &interpolation_domain) override
        {
            this->updateInterpolationDomain(interpolation_domain);
            fitCoefficientMatrix(data_matrix);
        }

        /// Evaluates the recursive Chebyshev basis terms at a scaled point.
        [[nodiscard]] vector_type evaluatePolynomialBase(const scalar_type scaled_eval_point) const override
        {
            vector_type polynomial_terms = vector_type::Zero(this->polynomial_degree_ + 1);
            polynomial_terms(0) = scalar_type(1);

            if (this->polynomial_degree_ > 0)
            {
                polynomial_terms(1) = scaled_eval_point;
            }

            for (int degree = 2; degree <= this->polynomial_degree_; ++degree)
            {
                polynomial_terms(degree) = scalar_type(2) * scaled_eval_point * polynomial_terms(degree - 1) -
                                           polynomial_terms(degree - 2);
            }

            return polynomial_terms;
        }

        /// Evaluates the fitted interpolant at the requested domain point.
        [[nodiscard]] vector_type evaluate(const scalar_type eval_point, const bool scale_value = true) const override
        {
            if (!this->checkEvalPointValidity(eval_point))
            {
                throw std::invalid_argument("Evaluation point is not within the interpolation domain.");
            }

            if (!this->hasCoefficients())
            {
                throw std::runtime_error("Polynomial coefficient matrix is empty.");
            }

            scalar_type scaled_eval_point = eval_point;
            if (scale_value)
            {
                // Evaluation uses the same affine scaling as fitting so the recursive basis sees
                // a point inside the canonical Chebyshev interval.
                scaled_eval_point = (scalar_type(2) * eval_point - (this->time_upper_bound_ + this->time_lower_bound_)) /
                                    (this->time_upper_bound_ - this->time_lower_bound_);
            }

            const auto polynomial_terms = evaluatePolynomialBase(scaled_eval_point);
            return this->coefficient_matrix_ * polynomial_terms;
        }
    };

    /// Chebyshev interpolator specialized for unit quaternion trajectories.
    template <TQuatLike Quaternion = Eigen::Vector4d>
    class CChebyshevQuaternionInterpolator
        : public CChebyshevInterpolator<quaternion_scalar_t<Quaternion>>
    {
        // TYPEDEFS
      public:
        using scalar_type = quaternion_scalar_t<Quaternion>;
        using base_type = CChebyshevInterpolator<scalar_type>;
        using typename base_type::matrix_type;
        using typename base_type::vector_type;

        // CONSTRUCTORS
        CChebyshevQuaternionInterpolator() = default;

        /// Constructs a quaternion interpolator without fitting it.
        /// @param polynomial_degree Degree of the fitted Chebyshev polynomial.
        /// @param interpolation_domain Domain associated with the quaternion samples.
        /// @param enable_sign_switch_patching When `true`, every fit removes consecutive sign flips before solving.
        CChebyshevQuaternionInterpolator(const int polynomial_degree,
                                         const vector_type &interpolation_domain,
                                         const bool enable_sign_switch_patching = true)
            : base_type(polynomial_degree, interpolation_domain, 4, true),
              enable_sign_switch_patching_(enable_sign_switch_patching)
        {
        }

        /// Constructs and fits a quaternion interpolator.
        /// @param polynomial_degree Degree of the fitted Chebyshev polynomial.
        /// @param interpolation_domain Domain associated with the quaternion samples.
        /// @param data_matrix Quaternion samples stored column-wise in scalar-first order.
        /// @param enable_sign_switch_patching When `true`, consecutive sign flips are removed before fitting.
        CChebyshevQuaternionInterpolator(const int polynomial_degree,
                                         const vector_type &interpolation_domain,
                                         const matrix_type &data_matrix,
                                         const bool enable_sign_switch_patching = true)
            : CChebyshevQuaternionInterpolator(polynomial_degree,
                                               interpolation_domain,
                                               enable_sign_switch_patching)
        {
            fitCoefficientMatrix(data_matrix);
        }

        /// Backwards-compatible constructor that keeps the historical parameter ordering.
        CChebyshevQuaternionInterpolator(const int polynomial_degree,
                                         const vector_type &interpolation_domain,
                                         const bool enable_sign_switch_patching,
                                         const matrix_type &data_matrix)
            : CChebyshevQuaternionInterpolator(polynomial_degree,
                                               interpolation_domain,
                                               data_matrix,
                                               enable_sign_switch_patching)
        {
        }

        // METHODS

        /// Fits the quaternion coefficient matrix for the currently configured domain.
        void fitCoefficientMatrix(const matrix_type &data_matrix)
        {
            validateQuaternionDataMatrix(data_matrix);

            if (!enable_sign_switch_patching_)
            {
                fitPreparedCoefficientMatrix(data_matrix);
                return;
            }

            matrix_type patched_data = data_matrix;
            this->patchDiscontinuities(patched_data);
            fitPreparedCoefficientMatrix(patched_data);
        }

        /// Fits the quaternion coefficient matrix after updating the interpolation domain.
        void fitCoefficientMatrix(const matrix_type &data_matrix, const vector_type &interpolation_domain) override
        {
            this->updateInterpolationDomain(interpolation_domain);
            fitCoefficientMatrix(data_matrix);
        }

        /// Evaluates the interpolant and returns a normalized quaternion value.
        [[nodiscard]] auto evaluateQuaternion(const scalar_type eval_point, const bool scale_value = true) const
            -> std::remove_cvref_t<Quaternion>
        {
            // Polynomial interpolation does not preserve unit norm, so renormalize before handing
            // the coefficients back to quaternion-aware code.
            const auto interpolated = base_type::evaluate(eval_point, scale_value).normalized();
            Eigen::Matrix<scalar_type, 4, 1> coefficients = interpolated.template head<4>();
            return makeQuatScalarFirst<Quaternion>(coefficients);
        }

        /// Evaluates the interpolant and returns the corresponding rotation matrix.
        [[nodiscard]] auto evaluateDcm(const scalar_type eval_point, const bool scale_value = true) const
            -> Eigen::Matrix<scalar_type, 3, 3>
        {
            return quaternionToDcm(evaluateQuaternion(eval_point, scale_value));
        }

        /// Checks interpolation accuracy against reference quaternion samples.
        /// @param data_points Ground-truth quaternion samples.
        /// @param evaluation_points Domain locations associated with `data_points`.
        /// @param absolute_tolerance Maximum allowed quaternion distance.
        /// @param scale_value Whether to scale `evaluation_points` before evaluation.
        /// @return `true` if all tested quaternion distances are below `absolute_tolerance`.
        [[nodiscard]] bool checkInterpolationError(const std::vector<vector_type> &data_points,
                                                   const vector_type &evaluation_points,
                                                   const scalar_type absolute_tolerance,
                                                   const bool scale_value) const
        {
            if (!this->hasCoefficients())
            {
                throw std::runtime_error("Polynomial coefficient matrix is empty.");
            }

            if (static_cast<int>(data_points.size()) != evaluation_points.size())
            {
                throw std::invalid_argument("Data points and evaluation points must have the same length.");
            }

            vector_type errors = vector_type::Zero(evaluation_points.size());
            for (Eigen::Index index = 0; index < evaluation_points.size(); ++index)
            {
                const auto interpolated = evaluateQuaternion(evaluation_points(index), scale_value);
                Eigen::Matrix<scalar_type, 4, 1> ground_truth = data_points[static_cast<std::size_t>(index)].template head<4>();
                errors(index) = quaternionDistance(interpolated, makeQuatScalarFirst<Quaternion>(ground_truth));
            }

            return errors.maxCoeff() < absolute_tolerance;
        }

      private:
        // PRIVATE METHODS / DATA MEMBERS
        /// Validates the quaternion sample matrix before preprocessing or fitting.
        void validateQuaternionDataMatrix(const matrix_type &data_matrix) const
        {
            if (data_matrix.rows() != 4)
            {
                throw std::invalid_argument("Quaternion interpolation expects a 4xN data matrix.");
            }
        }

        /// Converts a column-major quaternion data matrix into the vector list expected by the error checker.
        [[nodiscard]] static auto buildQuaternionDataPoints(const matrix_type &data_matrix)
            -> std::vector<vector_type>
        {
            std::vector<vector_type> data_points(static_cast<std::size_t>(data_matrix.cols()));
            for (Eigen::Index index = 0; index < data_matrix.cols(); ++index)
            {
                data_points[static_cast<std::size_t>(index)] = data_matrix.col(index);
            }

            return data_points;
        }

        /// Fits the prepared quaternion data and runs the quaternion-specific validation once.
        void fitPreparedCoefficientMatrix(const matrix_type &prepared_data)
        {
            const bool auto_check_enabled = this->auto_check_;

            // The base-class check uses Euclidean vector error, which is not the right metric for
            // unit quaternions. Disable it temporarily and run the quaternion-specific check below.
            this->setAutoCheck(false);
            try
            {
                base_type::fitCoefficientMatrix(prepared_data);
            }
            catch (...)
            {
                this->setAutoCheck(auto_check_enabled);
                throw;
            }

            this->setAutoCheck(auto_check_enabled);

            if (!auto_check_enabled)
            {
                return;
            }

            const auto data_points = buildQuaternionDataPoints(prepared_data);
            if (!checkInterpolationError(data_points, this->interpolation_domain_, scalar_type(1e-5), true))
            {
                throw std::runtime_error("Quaternion interpolation polynomial fitting failed.");
            }
        }

        bool enable_sign_switch_patching_{true};
    };
} // namespace mathcore
