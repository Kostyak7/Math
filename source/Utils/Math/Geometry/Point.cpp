#include "Point.h"

#include <cmath>

math::geom::Point2D::Point2D(value_type x, value_type y) 
    : x(x), y(y)
{
}

math::geom::Point2D::value_type math::geom::Point2D::distance_to(const Point2D& other) {
    value_type dx = x - other.x;
    value_type dy = y - other.y;
    return std::sqrt(dx * dx + dy * dy);
}

math::geom::Point2D math::geom::operator+(const Point2D& p1, const Point2D& p2) {
    Point2D res = p1;
    res.x += p2.x;
    res.y += p2.y;
    return res;
}

math::geom::Point2D math::geom::operator-(const Point2D& p1, const Point2D& p2) {
    Point2D res = p1;
    res.x -= p2.x;
    res.y -= p2.y;
    return res;
}

math::geom::Point2D math::geom::operator*(const Point2D& point, Point2D::value_type scalar) {
    Point2D res = point;
    res.x *= scalar;
    res.y *= scalar;
    return res;
}

math::geom::Point2D math::geom::operator*(Point2D::value_type scalar, const Point2D& point) {
    return point * scalar;
}

math::geom::Point2D math::geom::operator/(const Point2D& point, Point2D::value_type scalar) {
    if (scalar == 0) throw std::invalid_argument("Division by zero");
    return point * (1. / scalar);
}

math::geom::Point3D::Point3D(value_type x, value_type y, value_type z) 
    : x(x), y(y), z(z) 
{
}

math::geom::Point3D::value_type math::geom::Point3D::distance_to(const Point3D& other) {
    value_type dx = x - other.x;
    value_type dy = y - other.y;
    value_type dz = z - other.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

math::geom::Point3D math::geom::operator+(const Point3D& p1, const Point3D& p2) {
    Point3D res = p1;
    res.x += p2.x;
    res.y += p2.y;
    res.z += p2.z;
    return res;
}

math::geom::Point3D math::geom::operator-(const Point3D& p1, const Point3D& p2) {
    Point3D res = p1;
    res.x -= p2.x;
    res.y -= p2.y;
    res.z -= p2.z;
    return res;
}

math::geom::Point3D math::geom::operator*(const Point3D& point, Point3D::value_type scalar) {
    Point3D res = point;
    res.x *= scalar;
    res.y *= scalar;
    res.z *= scalar;
    return res;
}

math::geom::Point3D math::geom::operator*(Point3D::value_type scalar, const Point3D& point) {
    return point * scalar;
}

math::geom::Point3D math::geom::operator/(const Point3D& point, Point3D::value_type scalar) {
    if (scalar == 0) throw std::invalid_argument("Division by zero");
    return point * (1. / scalar);
}

math::geom::Point4D::Point4D(value_type x, value_type y, value_type z, value_type w) 
    : x(x), y(y), z(z), w(w) 
{
}

math::geom::Point4D::value_type math::geom::Point4D::distance_to(const Point4D& other) {
    value_type dx = x - other.x;
    value_type dy = y - other.y;
    value_type dz = z - other.z;
    value_type dw = w - other.w;
    return std::sqrt(dx * dx + dy * dy + dz * dz + dw * dw);
}

math::geom::Point4D math::geom::operator+(const Point4D& p1, const Point4D& p2) {
    Point4D res = p1;
    res.x += p2.x;
    res.y += p2.y;
    res.z += p2.z;
    res.w += p2.w;
    return res;
}

math::geom::Point4D math::geom::operator-(const Point4D& p1, const Point4D& p2) {
    Point4D res = p1;
    res.x -= p2.x;
    res.y -= p2.y;
    res.z -= p2.z;
    res.w -= p2.w;
    return res;
}

math::geom::Point4D math::geom::operator*(const Point4D& point, Point4D::value_type scalar) {
    Point4D res = point;
    res.x *= scalar;
    res.y *= scalar;
    res.z *= scalar;
    res.w *= scalar;
    return res;
}

math::geom::Point4D math::geom::operator*(Point4D::value_type scalar, const Point4D& point) {
    return point * scalar;
}

math::geom::Point4D math::geom::operator/(const Point4D& point, Point4D::value_type scalar) {
    if (scalar == 0) throw std::invalid_argument("Division by zero");
    return point * (1. / scalar);
}