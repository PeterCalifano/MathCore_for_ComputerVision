#pragma once

/// @file
/// @brief Traits and concepts used to adapt quaternion-like types.

#include <array>
#include <concepts>
#include <type_traits>

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace mathcore
{
    /// Tag type for quaternions stored as `[w, x, y, z]`.
    struct SScalarFirst
    {
    };

    /// Tag type for quaternions stored as `[x, y, z, w]`.
    struct SScalarLast
    {
    };

    /// Concept for supported quaternion storage conventions.
    template <class T>
    concept TQuatConvention = std::same_as<T, SScalarFirst> || std::same_as<T, SScalarLast>;

    /**
     * @brief Primary customization point for integrating quaternion types with MathCore. Specialize this struct for a custom quaternion type to enable compatibility with MathCore's quaternion utilities and interpolators.
     *
     * @tparam T
     * @tparam typename
     */
    template <typename T, typename = void>
    struct SQuatTraits;

    /**
     * @brief Specialization of `SQuatTraits` for `std::array<Scalar, 4>` with scalar-first convention.
     *
     * @tparam Scalar
     */
    template <typename Scalar>
    struct SQuatTraits<std::array<Scalar, 4>>
    {
        using scalar_type = Scalar;
        using convention = SScalarFirst;

        static auto toScalarFirst(const std::array<Scalar, 4> &value) -> Eigen::Matrix<scalar_type, 4, 1>
        {
            Eigen::Matrix<scalar_type, 4, 1> coefficients;
            coefficients << value[0], value[1], value[2], value[3];
            return coefficients;
        }

        static auto fromScalarFirst(const Eigen::Matrix<scalar_type, 4, 1> &coefficients) -> std::array<Scalar, 4>
        {
            return {coefficients(0), coefficients(1), coefficients(2), coefficients(3)};
        }
    };

    /**
     * @brief Specialization of `SQuatTraits` for `Eigen::Quaternion` types with scalar-first convention.
     *
     * @tparam Scalar
     * @tparam Options
     */
    template <typename Scalar, int Options>
    struct SQuatTraits<Eigen::Quaternion<Scalar, Options>>
    {
        using scalar_type = Scalar;
        using convention = SScalarFirst;

        static auto toScalarFirst(const Eigen::Quaternion<Scalar, Options> &value) -> Eigen::Matrix<scalar_type, 4, 1>
        {
            Eigen::Matrix<scalar_type, 4, 1> coefficients;
            coefficients << value.w(), value.x(), value.y(), value.z();
            return coefficients;
        }

        static auto fromScalarFirst(const Eigen::Matrix<scalar_type, 4, 1> &coefficients) -> Eigen::Quaternion<Scalar, Options>
        {
            return {coefficients(0), coefficients(1), coefficients(2), coefficients(3)};
        }
    };

    /**
     * @brief Specialization of `SQuatTraits` for 4x1 Eigen vectors with scalar-first convention.
     * @details This allows raw Eigen vectors to be used as quaternions in MathCore, with the expectation that they are ordered as `[w, x, y, z]`.
     *
     * @tparam Scalar
     * @tparam Options
     * @tparam MaxRows
     * @tparam MaxCols
     */
    template <typename Scalar, int Options, int MaxRows, int MaxCols>
    struct SQuatTraits<Eigen::Matrix<Scalar, 4, 1, Options, MaxRows, MaxCols>>
    {
        using scalar_type = Scalar;
        using convention = SScalarFirst;

        static auto toScalarFirst(const Eigen::Matrix<Scalar, 4, 1, Options, MaxRows, MaxCols> &value)
            -> Eigen::Matrix<scalar_type, 4, 1>
        {
            return value;
        }

        static auto fromScalarFirst(const Eigen::Matrix<scalar_type, 4, 1> &coefficients)
            -> Eigen::Matrix<Scalar, 4, 1, Options, MaxRows, MaxCols>
        {
            return coefficients;
        }
    };

    /**
     * @brief Alias for the `SQuatTraits` specialization associated with a quaternion-like type `T`.
     * @details Used internally to simplify access to the traits within MathCore's quaternion utilities and interpolators.
     *
     * @tparam T
     */
    template <typename T>
    using SQuatTraits_t = SQuatTraits<std::remove_cvref_t<T>>;

    /// Concept for quaternion types adapted through `SQuatTraits`

    /**
     * @brief Concept for quaternion-like types that can be adapted through `SQuatTraits`
     * @note To satisfy this concept, a type `T` must have a corresponding specialization of `SQuatTraits` that defines a scalar type, a storage convention, and provides methods for converting to and from scalar-first coefficient vectors. This allows MathCore's quaternion utilities and interpolators to operate on a wide range of quaternion representations as long as they are properly adapted through `SQuatTraits`.
     * 
     * @tparam T 
     */
    template <typename T>
    concept TQuatLike = requires(const std::remove_cvref_t<T> &value,
                                 const Eigen::Matrix<typename SQuatTraits_t<T>::scalar_type, 4, 1> &coefficients) {
        typename SQuatTraits_t<T>::scalar_type;
        typename SQuatTraits_t<T>::convention;
        requires TQuatConvention<typename SQuatTraits_t<T>::convention>;
        { SQuatTraits_t<T>::toScalarFirst(value) } -> std::same_as<Eigen::Matrix<typename SQuatTraits_t<T>::scalar_type, 4, 1>>;
        { SQuatTraits_t<T>::fromScalarFirst(coefficients) } -> std::same_as<std::remove_cvref_t<T>>;
    };

    /**
     * @brief Scalar type associated with a quaternion-like type `Quaternion`, as defined by its `SQuatTraits` specialization. This is the type of the individual coefficients of the quaternion.
     * 
     * @tparam Quaternion 
     */
    template <TQuatLike Quaternion>
    using quaternion_scalar_t = typename SQuatTraits_t<Quaternion>::scalar_type;

    /**
     * @brief Storage convention associated with a quaternion-like type `Quaternion`, as defined by its `SQuatTraits` specialization. This indicates whether the quaternion is stored in scalar-first (`[w, x, y, z]`) or scalar-last (`[x, y, z, w]`) order.
     * 
     * @tparam Quaternion 
     */
    template <TQuatLike Quaternion>
    using quaternion_convention_t = typename SQuatTraits_t<Quaternion>::convention;

    /**
     * @brief Type alias for a 4x1 Eigen vector containing the coefficients of a quaternion-like type `Quaternion` in scalar-first order.
     * 
     * @tparam Quaternion 
     */
    template <TQuatLike Quaternion>
    using quaternion_coefficients_t = Eigen::Matrix<quaternion_scalar_t<Quaternion>, 4, 1>;

    /**
     * @brief Extracts the coefficients of a quaternion-like type `Quaternion` in scalar-first order.
     * 
     * @tparam Quaternion 
     * @param value 
     * @return quaternion_coefficients_t<Quaternion> 
     */
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto getCoefficientsScalarFirst(const Quaternion &value)
        -> quaternion_coefficients_t<Quaternion>
    {
        return SQuatTraits_t<Quaternion>::toScalarFirst(value);
    }

    /**
     * @brief Reconstructs a quaternion-like object from scalar-first coefficients.
     * 
     * @tparam Quaternion 
     * @param coefficients 
     * @return std::remove_cvref_t<Quaternion> 
     */
    template <TQuatLike Quaternion>
    [[nodiscard]] inline auto makeQuatScalarFirst(const quaternion_coefficients_t<Quaternion> &coefficients)
        -> std::remove_cvref_t<Quaternion>
    {
        return SQuatTraits_t<Quaternion>::fromScalarFirst(coefficients);
    }
} // namespace mathcore
