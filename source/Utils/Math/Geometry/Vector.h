#pragma once

#include "Point.h"

namespace math::geom {

    class Vector2D {
    public:
        using value_type = Point2D::value_type;

        value_type x, y;

        Vector2D(value_type x = 0.0, value_type y = 0.0);
        Vector2D(const Point2D& to);
        Vector2D(const Point2D& from, const Point2D& to);

        const value_type& operator[](size_t index) const;
        value_type& operator[](size_t index);

        Vector2D normalized() const;
        value_type norm() const;
        value_type length() const;
        value_type dot(const Vector2D& other) const;
        value_type pseudodot(const Vector2D& other) const;
    };

    Vector2D operator+(const Vector2D& v1, const Vector2D& v2);
    Vector2D operator-(const Vector2D& v1, const Vector2D& v2);
    Vector2D operator*(const Vector2D& vec, Vector2D::value_type scalar);
    Vector2D operator*(Vector2D::value_type scalar, const Vector2D& vec);
    Vector2D operator/(const Vector2D& vec, Vector2D::value_type scalar);

    class Vector3D {
    public:
        using value_type = Point3D::value_type;

        value_type x, y, z;

        Vector3D(value_type x = 0.0, value_type y = 0.0, value_type z = 0.0);
        Vector3D(const Point3D& to);
        Vector3D(const Point3D& from, const Point3D& to);

        const value_type& operator[](size_t index) const;
        value_type& operator[](size_t index);

        Vector3D normalized() const;
        value_type norm() const;
        value_type length() const;
        value_type dot(const Vector3D& other) const;
        value_type cross(const Vector3D& other) const;
    };

    Vector3D operator+(const Vector3D& v1, const Vector3D& v2);
    Vector3D operator-(const Vector3D& v1, const Vector3D& v2);
    Vector3D operator*(const Vector3D& vec, Vector3D::value_type scalar);
    Vector3D operator*(Vector3D::value_type scalar, const Vector3D& vec);
    Vector3D operator/(const Vector3D& vec, Vector3D::value_type scalar);

    class Vector4D {
    public:
        using value_type = Point4D::value_type;

        value_type x, y;

        Vector4D(value_type x = 0.0, value_type y = 0.0, value_type z = 0.0, value_type w = 0.0);
        Vector4D(const Point4D& to);
        Vector4D(const Point4D& from, const Point4D& to);

        const value_type& operator[](size_t index) const;
        value_type& operator[](size_t index);

        Vector4D normalized() const;
        value_type norm() const;
        value_type length() const;
        value_type dot(const Vector4D& other) const;
    };

    Vector4D operator+(const Vector4D& v1, const Vector4D& v2);
    Vector4D operator-(const Vector4D& v1, const Vector4D& v2);
    Vector4D operator*(const Vector4D& vec, Vector4D::value_type scalar);
    Vector4D operator*(Vector4D::value_type scalar, const Vector4D& vec);
    Vector4D operator/(const Vector4D& vec, Vector4D::value_type scalar);

} // namespace math::geom