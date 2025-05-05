#pragma once

namespace math::geom {

    struct Point2D {
        using value_type = double;

        value_type x, y;

        Point2D(value_type x = 0.0, value_type y = 0.0);
        value_type distance_to(const Point2D& other);
    };

    Point2D operator+(const Point2D& p1, const Point2D& p2);
    Point2D operator-(const Point2D& p1, const Point2D& p2);
    Point2D operator*(const Point2D& point, Point2D::value_type scalar);
    Point2D operator*(Point2D::value_type scalar, const Point2D& point);
    Point2D operator/(const Point2D& point, Point2D::value_type scalar);

    struct Point3D {
        using value_type = double;

        value_type x, y, z;

        Point3D(value_type x = 0.0, value_type y = 0.0, value_type z = 0.0);
        value_type distance_to(const Point3D& other);
    };

    Point3D operator+(const Point3D& p1, const Point3D& p2);
    Point3D operator-(const Point3D& p1, const Point3D& p2);
    Point3D operator*(const Point3D& point, Point3D::value_type scalar);
    Point3D operator*(Point3D::value_type scalar, const Point3D& point);
    Point3D operator/(const Point3D& point, Point3D::value_type scalar);

    struct Point4D {
        using value_type = double;

        value_type x, y, z, w;

        Point4D(value_type x = 0.0, value_type y = 0.0, value_type z = 0.0, value_type w = 0.0);
        value_type distance_to(const Point4D& other);
    };

    Point4D operator+(const Point4D& p1, const Point4D& p2);
    Point4D operator-(const Point4D& p1, const Point4D& p2);
    Point4D operator*(const Point4D& point, Point4D::value_type scalar);
    Point4D operator*(Point4D::value_type scalar, const Point4D& point);
    Point4D operator/(const Point4D& point, Point4D::value_type scalar);


} // namespace math::geom
