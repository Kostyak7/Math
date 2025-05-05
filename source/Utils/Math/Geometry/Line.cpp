#include "Line.h"

#include <Utils/Math/Common.h>

#include <cmath>
#include <stdexcept>

math::geom::Line2D::Line2D(const Vector2D& vec)
    : p1(0., 0.), p2(vec.x, vec.y) 
{
}

math::geom::Line2D::Line2D(const Point2D& point)
    : p1(0., 0.), p2(point.x, point.y) 
{
}

math::geom::Line2D::Line2D(const Point2D& p1, const Point2D& p2)
    : p1(p1), p2(p2) 
{
}

math::geom::Point2D::value_type math::geom::Line2D::length() const {
    return std::hypot(p2.x - p1.x, p2.y - p1.y);
}

bool math::geom::Line2D::contains(const Point2D& point) const {
    Vector2D v1(p1, p2);
    Vector2D v2(p1, point);
    return fabs(v1.pseudodot(v2)) < std::numeric_limits<Point2D::value_type>::epsilon();
}

bool math::geom::Line2D::parallel(const Line2D& other) const {
    Vector2D v1(p1, p2);
    Vector2D v2(other.p1, other.p2);
    return fabs(v1.pseudodot(v2)) < std::numeric_limits<Point2D::value_type>::epsilon();
}

bool math::geom::Line2D::collinear(const Line2D& other) const {
    return parallel(other) && contains(other.p1);
}

bool math::geom::Line2D::perpendicular(const Line2D& other) const {
    Vector2D v1(p1, p2);
    Vector2D v2(other.p1, other.p2);

    return fabs(v1.dot(v2)) < std::numeric_limits<Point2D::value_type>::epsilon();
}

math::geom::Line3D::Line3D(const Vector3D& vec)
    : p1(0., 0., 0.), p2(vec.x, vec.y, vec.z) 
{
}

math::geom::Line3D::Line3D(const Point3D& point)
    : p1(0., 0., 0.), p2(point.x, point.y, point.z)
{
}

math::geom::Line3D::Line3D(const Point3D& p1, const Point3D& p2)
    : p1(p1), p2(p2) 
{
}

math::geom::Point3D::value_type math::geom::Line3D::length() const {
    return std::sqrt(
        std::pow(p2.x - p1.x, 2) +
        std::pow(p2.y - p1.y, 2) +
        std::pow(p2.z - p1.z, 2));
}

bool math::geom::Line3D::contains(const Point3D& point) const {
    Vector3D line_vec(p1, p2);
    Vector3D point_vec(p1, point);

    Vector3D cross = line_vec.cross(point_vec);
    return fabs(cross.x) < std::numeric_limits<Point3D::value_type>::epsilon() &&
        fabs(cross.y) < std::numeric_limits<Point3D::value_type>::epsilon() &&
        fabs(cross.z) < std::numeric_limits<Point3D::value_type>::epsilon();
}

bool math::geom::Line3D::parallel(const Line3D& other) const {
    Vector3D v1(p1, p2);
    Vector3D v2(other.p1, other.p2);

    Vector3D cross = v1.cross(v2);
    return fabs(cross.x) < std::numeric_limits<Point3D::value_type>::epsilon() &&
        fabs(cross.y) < std::numeric_limits<Point3D::value_type>::epsilon() &&
        fabs(cross.z) < std::numeric_limits<Point3D::value_type>::epsilon();
}

bool math::geom::Line3D::collinear(const Line3D& other) const {
    return parallel(other) && contains(other.p1);
}

bool math::geom::Line3D::perpendicular(const Line3D& other) const {
    Vector3D v1(p1, p2);
    Vector3D v2(other.p1, other.p2);

    return fabs(v1.dot(v2)) < std::numeric_limits<Point3D::value_type>::epsilon();
}