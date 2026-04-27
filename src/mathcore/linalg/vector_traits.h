#pragma once

/// @file
/// @brief Traits and concepts used to adapt 3D vector-like types.

#include <concepts>
#include <cstddef>
#include <type_traits>

#include <Eigen/Dense>

namespace mathcore
{
    /// Primary customization point for integrating external vector types with MathCore.
    template <typename T, typename = void>
    struct vector_traits;

    /// Convenience alias for the traits specialization associated with `T`.
    template <typename T>
    using vector_traits_t = vector_traits<std::remove_cvref_t<T>>;

    /// Concept for fixed-size vector types that expose MathCore vector traits.
    template <typename T>
    concept vector_like = requires(const std::remove_cvref_t<T> &value,
                                   const Eigen::Matrix<typename vector_traits_t<T>::scalar_type,
                                                       static_cast<Eigen::Index>(vector_traits_t<T>::static_size),
                                                       1> &coefficients,
                                   const std::size_t index) {
        typename vector_traits_t<T>::scalar_type;
        requires vector_traits_t<T>::is_specialized;
        { vector_traits_t<T>::get(value, index) } -> std::convertible_to<typename vector_traits_t<T>::scalar_type>;
        { vector_traits_t<T>::make(coefficients) } -> std::same_as<std::remove_cvref_t<T>>;
        { vector_traits_t<T>::zero() } -> std::same_as<std::remove_cvref_t<T>>;
    };

    /// Scalar type associated with a vector-like type.
    template <vector_like T>
    using vector_scalar_t = typename vector_traits_t<T>::scalar_type;

    /// Compile-time vector size associated with a vector-like type.
    template <vector_like T>
    inline constexpr std::size_t vector_size_v = vector_traits_t<T>::static_size;

    /// Concept for vector-like types with exactly three components.
    template <typename T>
    concept vector3_like = vector_like<T> && (vector_size_v<T> == 3);

    /// MathCore adapter for Eigen 3x1 column vectors.
    template <typename Scalar, int Options, int MaxRows, int MaxCols>
    struct vector_traits<Eigen::Matrix<Scalar, 3, 1, Options, MaxRows, MaxCols>>
    {
        // TYPEDEFS
        using scalar_type = Scalar;
        using vector_type = Eigen::Matrix<Scalar, 3, 1, Options, MaxRows, MaxCols>;

        // CONSTANTS
        static constexpr bool is_specialized = true;
        static constexpr std::size_t static_size = 3;

        // METHODS
        static scalar_type get(const vector_type &value, const std::size_t index)
        {
            return value(static_cast<Eigen::Index>(index));
        }

        static vector_type make(const Eigen::Matrix<Scalar, 3, 1> &coefficients)
        {
            return coefficients;
        }

        static vector_type zero()
        {
            return vector_type::Zero();
        }
    };

    /// Converts an arbitrary supported 3D vector type into an Eigen column vector.
    /// @param value Source vector expressed through `vector_traits`.
    /// @return A 3x1 Eigen vector with the same coefficients as `value`.
    template <vector3_like Vector>
    [[nodiscard]] inline auto to_eigen_vector3(const Vector &value)
        -> Eigen::Matrix<vector_scalar_t<Vector>, 3, 1>
    {
        using Scalar = vector_scalar_t<Vector>;

        Eigen::Matrix<Scalar, 3, 1> result;
        for (std::size_t index = 0; index < 3; ++index)
        {
            result(static_cast<Eigen::Index>(index)) = vector_traits_t<Vector>::get(value, index);
        }

        return result;
    }

} // namespace mathcore
