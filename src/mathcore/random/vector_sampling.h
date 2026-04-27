#pragma once
#include <mathcore/linalg/vector_traits.h>
#include <random>

namespace mathcore
{
    /**
     * @brief Draws a 3D Gaussian perturbation around a mean vector.
     *
     * @tparam Vector
     * @tparam Generator
     * @param standard_deviation
     * @param mean
     * @param seed
     * @return std::remove_cvref_t<Vector>
     */
    template <vector3_like Vector = Eigen::Vector3d, typename Generator = std::mt19937>
        requires std::floating_point<vector_scalar_t<Vector>>
    [[nodiscard]] inline auto scatter_vector3(const vector_scalar_t<Vector> standard_deviation = vector_scalar_t<Vector>(0),
                                              const Vector &mean = vector_traits_t<Vector>::zero(),
                                              const typename Generator::result_type seed = 55) -> std::remove_cvref_t<Vector>
    {
        if (standard_deviation <= vector_scalar_t<Vector>(0))
        {
            // Return mean if standard deviation is zero or negative
            return mean;
        }

        // Generator and distribution setup
        Generator generator(seed);
        std::normal_distribution<vector_scalar_t<Vector>> distribution(vector_scalar_t<Vector>(0), standard_deviation);

        // Allocate temporary Eigen vector to hold the sampled coefficients before converting back to the target type
        Eigen::Matrix<vector_scalar_t<Vector>, 3, 1> sampled;
        sampled << distribution(generator), distribution(generator), distribution(generator);

        return vector_traits_t<Vector>::make(sampled + to_eigen_vector3(mean));
    }
} // namespace mathcore
