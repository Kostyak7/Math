#include <Math/Geometry/Intersections.h>

#include <stdexcept>

bool math::geom::intersect(const Line3D& line, const Plane& plane, Point3D* out_point) {
    // ...
    return {};
}

std::optional<math::geom::Point2D> math::geom::intersection_point(const Line2D& l1, const Line2D& l2) {
    // ...
    return {};
}

std::optional<math::geom::Point3D> math::geom::intersection_point(const Line3D& l1, const Line3D& l2) {
    // ...
    return {};
}

std::optional<math::geom::Point3D> math::geom::intersection_point(const Line3D& line, const Plane& plane) {
    // ...
    return {};
}

std::optional<math::geom::Point2D> math::geom::intersection_point(const Vector2D& v1, const Vector2D& v2) {
    // ...
    return {};
}

std::optional<math::geom::Point3D> math::geom::intersection_point(const Vector3D& v1, const Vector3D& v2) {
    // ...
    return {};
}

std::optional<math::geom::Point3D> math::geom::intersection_point(const Vector3D& vec, const Plane& plane) {
    // ...
    return {};
}

std::optional<math::geom::Vector3D> math::geom::intersection_vector(const Plane& pln1, const Plane& pln2) {
    // ...
    return {};
}

bool math::geom::is_point_in_triangle(const Point2D& p, const Point2D& a, const Point2D& b, const Point2D& c) {
    double area = triangle_area(a, b, c);
    double area1 = triangle_area(p, b, c);
    double area2 = triangle_area(a, p, c);
    double area3 = triangle_area(a, b, p);
    return is_near(area1 + area2 + area3, area);
}

bool math::geom::is_point_in_tetrahedron(const Point3D& p, const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
    double v_total = std::abs(tetrahedron_volume(a, b, c, d));
    double v1 = std::abs(tetrahedron_volume(p, b, c, d));
    double v2 = std::abs(tetrahedron_volume(a, p, c, d));
    double v3 = std::abs(tetrahedron_volume(a, b, p, d));
    double v4 = std::abs(tetrahedron_volume(a, b, c, p));

    return is_near(v1 + v2 + v3 + v4, v_total);
}
