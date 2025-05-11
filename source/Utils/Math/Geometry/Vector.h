#pragma once

#include "Point.h"

namespace math::geom {

    template <size_t N>
    struct Vector : public std::array<double, N> {
        using value_type = double;
        using base_type = std::array<value_type, N>;
        using base_type::base_type;
        Vector(const Point<N, value_type>& to);
        Vector(const Point<N, value_type>& from, const Point<N, value_type>& to);
        Vector(std::initializer_list<value_type> init);

        Vector& operator*=(value_type scalar);
        Vector& operator/=(value_type scalar);

        Vector& operator+=(const Vector& other);
        Vector& operator-=(const Vector& other);

        bool is_zero() const;

        value_type x() const;
        value_type y() const;
        value_type z() const;
        value_type w() const;

        value_type dot(const Vector& other) const;
        value_type norm() const;
        value_type length() const;

        Vector& normalize();
    };

    template <size_t N>
    Vector<N> operator*(const Vector<N>& vec, typename Vector<N>::value_type scalar);
    template <size_t N>
    Vector<N> operator*(typename Vector<N>::value_type scalar, const Vector<N>& vec);
    template <size_t N>
    Vector<N> operator/(const Vector<N>& vec, typename Vector<N>::value_type scalar);

    template <size_t N>
    Vector<N> operator+(const Vector<N>& v1, const Vector<N>& v2);
    template <size_t N>
    Vector<N> operator-(const Vector<N>& v1, const Vector<N>& v2);

    template <size_t N>
    Vector<N> operator-(const Vector<N>& vec);

    template <size_t N>
    Vector<N> normalized(const Vector<N>& vector);

    template <size_t N>
    Vector<N> zero_vector();
    template <size_t N>
    Vector<N> unit_vector(size_t i, typename Vector<N>::value_type value = 1.0);

    using Vector2D = Vector<2>;
    using Vector3D = Vector<3>;
    using Vector4D = Vector<4>;

    inline Vector2D::value_type pseudodot(const Vector2D& v1, const Vector2D& v2);

    inline Vector3D cross_product(const Vector3D& v1, const Vector3D& v2);

} // namespace math::geom

#include "Vector.hpp"
