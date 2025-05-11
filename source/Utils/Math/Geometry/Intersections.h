#pragma once 

#include "Utils.h"

#include <optional>

namespace math::geom {

    bool intersect(const Line3D& line, const Plane& plane, Point3D* out_point = nullptr);

    std::optional<Point2D> intersection_point(const Line2D& l1, const Line2D& l2);
    std::optional<Point3D> intersection_point(const Line3D& l1, const Line3D& l2);
    std::optional<Point3D> intersection_point(const Line3D& line, const Plane& plane);

    std::optional<Point2D> intersection_point(const Vector2D& v1, const Vector2D& v2);
    std::optional<Point3D> intersection_point(const Vector3D& v1, const Vector3D& v2);
    std::optional<Point3D> intersection_point(const Vector3D& vec, const Plane& plane);

    std::optional<Vector3D> intersection_vector(const Plane& plane, const Plane& plane);

    bool is_point_in_triangle(const Point2D& p, const Point2D& a, const Point2D& b, const Point2D& c);
    bool is_point_in_tetrahedron(const Point3D& p, const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d);

    // bool is_point_in_convex_polygon(...)

} // namespace math::geom
