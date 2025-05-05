#pragma once

#include "Vector.h"

namespace math::geom {

    class Line2D {
    public:
        Point2D p1, p2;

        Line2D(const Vector2D& vec);
        Line2D(const Point2D& point);
        Line2D(const Point2D& p1, const Point2D& p2);

        Point2D::value_type length() const;
        bool contains(const Point2D& point) const;
        bool parallel(const Line2D& other) const;
        bool collinear(const Line2D& other) const;
        bool perpendicular(const Line2D& other) const;
    };

    class Line3D {
    public:
        Point3D p1, p2;

        Line3D(const Vector3D& vec);
        Line3D(const Point3D& point);
        Line3D(const Point3D& p1, const Point3D& p2);

        Point3D::value_type length() const;
        bool contains(const Point3D& point) const;
        bool parallel(const Line3D& other) const;
        bool collinear(const Line3D& other) const;
        bool perpendicular(const Line3D& other) const;
    };

} // namespace math::geom