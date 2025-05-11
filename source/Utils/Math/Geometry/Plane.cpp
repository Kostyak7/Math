#include "Plane.h"

#include <Utils/Math/Common.h>

#include <cmath>
#include <stdexcept>

math::geom::Plane::Plane(const Point3D& p1, const Point3D& p2, const Point3D& p3) {
    Vector3D v1(p1, p2);
    Vector3D v2(p1, p3);
    Vector3D normal = normalized(cross_product(v1, v2));

    a = normal.x;
    b = normal.y;
    c = normal.z;
    d = -(a * p1.x + b * p1.y + c * p1.z);
}

math::geom::Vector3D math::geom::Plane::normal_vector() const {
    return Vector3D{ a, b, c };
}

math::geom::Plane::value_type math::geom::Plane::distance_to(const Point3D& point) const {
    value_type numerator = std::abs(a * point.x + b * point.y + c * point.z + d);
    value_type denominator = std::sqrt(a * a + b * b + c * c);
    return numerator / denominator;
}

bool math::geom::Plane::contains(const Point3D& point) const {
    return std::abs(a * point.x + b * point.y + c * point.z + d) < std::numeric_limits<Point2D::value_type>::epsilon();
}

bool math::geom::Plane::contains(const Vector3D& vec) const {
    return std::abs(a * vec.x + b * vec.y + c * vec.z) < std::numeric_limits<Point2D::value_type>::epsilon();
}
