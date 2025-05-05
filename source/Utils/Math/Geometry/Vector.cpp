#include "Vector.h"

#include <Utils/Math/Common.h>

#include <cmath>
#include <stdexcept>

math::geom::Vector2D::Vector2D(value_type x, value_type y)
    : x(x), y(y) 
{
}

math::geom::Vector2D::Vector2D(const Point2D& to)
    : x(to.x), y(to.y) 
{
}

math::geom::Vector2D::Vector2D(const Point2D& from, const Point2D& to)
    : x(to.x - from.x), y(to.y - from.y) 
{
}

const math::geom::Vector2D::value_type& math::geom::Vector2D::operator[](size_t index) const {
    switch (index) {
    case 0: return x;
    case 1: return y;
    default: throw std::out_of_range("Vector2D index out of range");
    }
}

math::geom::Vector2D::value_type& math::geom::Vector2D::operator[](size_t index) {
    switch (index) {
    case 0: return x;
    case 1: return y;
    default: throw std::out_of_range("Vector2D index out of range");
    }
}

math::geom::Vector2D math::geom::Vector2D::normalized() const {
    if (value_type len = norm(); dcmp(len) != 0)
        return Vector2D(x / len, y / len);
    return Vector2D();
}

math::geom::Vector2D::value_type math::geom::Vector2D::norm() const {
    return std::sqrt(x * x + y * y);
}

math::geom::Vector2D::value_type math::geom::Vector2D::length() const {
    return norm();
}

math::geom::Vector2D::value_type math::geom::Vector2D::dot(const Vector2D& other) const {
    return x * other.x + y * other.y;
}

math::geom::Vector2D::value_type math::geom::Vector2D::pseudodot(const Vector2D& other) const {
    return x * other.y - y * other.x;
}

math::geom::Vector2D math::geom::operator+(const Vector2D& v1, const Vector2D& v2) {
    return Vector2D(v1.x + v2.x, v1.y + v2.y);
}

math::geom::Vector2D math::geom::operator-(const Vector2D& v1, const Vector2D& v2) {
    return Vector2D(v1.x - v2.x, v1.y - v2.y);
}

math::geom::Vector2D math::geom::operator*(const Vector2D& vec, Vector2D::value_type scalar) {
    return Vector2D(vec.x * scalar, vec.y * scalar);
}

math::geom::Vector2D math::geom::operator*(Vector2D::value_type scalar, const Vector2D& vec) {
    return vec * scalar;
}

math::geom::Vector2D math::geom::operator/(const Vector2D& vec, Vector2D::value_type scalar) {
    if (scalar == 0) throw std::invalid_argument("Division by zero");
    return vec * (1. / scalar);
}

math::geom::Vector3D::Vector3D(value_type x, value_type y, value_type z)
    : x(x), y(y), z(z) 
{
}

math::geom::Vector3D::Vector3D(const Point3D& to)
    : x(to.x), y(to.y), z(to.z) 
{
}

math::geom::Vector3D::Vector3D(const Point3D& from, const Point3D& to)
    : x(to.x - from.x), y(to.y - from.y), z(to.z - from.z) 
{
}

const math::geom::Vector3D::value_type& math::geom::Vector3D::operator[](size_t index) const {
    switch (index) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: throw std::out_of_range("Vector3D index out of range");
    }
}

math::geom::Vector3D::value_type& math::geom::Vector3D::operator[](size_t index) {
    switch (index) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: throw std::out_of_range("Vector3D index out of range");
    }
}

math::geom::Vector3D math::geom::Vector3D::normalized() const {
    if (value_type len = norm(); dcmp(len) != 0)
        return Vector3D(x / len, y / len, z / len);
    return Vector3D();
}

math::geom::Vector3D::value_type math::geom::Vector3D::norm() const {
    return std::sqrt(x * x + y * y + z * z);
}

math::geom::Vector3D::value_type math::geom::Vector3D::length() const {
    return norm();
}

math::geom::Vector3D::value_type math::geom::Vector3D::dot(const Vector3D& other) const {
    return x * other.x + y * other.y + z * other.z;
}

math::geom::Vector3D math::geom::Vector3D::cross(const Vector3D& other) const {
    return Vector3D(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    );
}

math::geom::Vector3D math::geom::operator+(const Vector3D& v1, const Vector3D& v2) {
    return Vector3D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

math::geom::Vector3D math::geom::operator-(const Vector3D& v1, const Vector3D& v2) {
    return Vector3D(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

math::geom::Vector3D math::geom::operator*(const Vector3D& vec, Vector3D::value_type scalar) {
    return Vector3D(vec.x * scalar, vec.y * scalar, vec.z * scalar);
}

math::geom::Vector3D math::geom::operator*(Vector3D::value_type scalar, const Vector3D& vec) {
    return vec * scalar;
}

math::geom::Vector3D math::geom::operator/(const Vector3D& vec, Vector3D::value_type scalar) {
    if (scalar == 0) throw std::invalid_argument("Division by zero");
    return vec * (1. / scalar);
}

math::geom::Vector4D::Vector4D(value_type x, value_type y, value_type z, value_type w)
    : x(x), y(y), z(z), w(w) 
{
}

math::geom::Vector4D::Vector4D(const Point4D& to)
    : x(to.x), y(to.y), z(to.z), w(to.w) 
{
}

math::geom::Vector4D::Vector4D(const Point4D& from, const Point4D& to)
    : x(to.x - from.x), y(to.y - from.y), z(to.z - from.z), w(to.w - from.w) 
{
}

const math::geom::Vector4D::value_type& math::geom::Vector4D::operator[](size_t index) const {
    switch (index) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    case 3: return w;
    default: throw std::out_of_range("Vector4D index out of range");
    }
}

math::geom::Vector4D::value_type& math::geom::Vector4D::operator[](size_t index) {
    switch (index) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    case 3: return w;
    default: throw std::out_of_range("Vector4D index out of range");
    }
}

math::geom::Vector4D math::geom::Vector4D::normalized() const {
    if (value_type len = norm(); dcmp(len) != 0)
        return Vector4D(x / len, y / len, z / len, w / len);
    return Vector4D();
}

math::geom::Vector4D::value_type math::geom::Vector4D::norm() const {
    return std::sqrt(x * x + y * y + z * z + w * w);
}

math::geom::Vector4D::value_type math::geom::Vector4D::length() const {
    return norm();
}

math::geom::Vector4D::value_type math::geom::Vector4D::dot(const Vector4D& other) const {
    return x * other.x + y * other.y + z * other.z + w * other.w;
}

math::geom::Vector4D math::geom::operator+(const Vector4D& v1, const Vector4D& v2) {
    return Vector4D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w);
}

math::geom::Vector4D math::geom::operator-(const Vector4D& v1, const Vector4D& v2) {
    return Vector4D(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w);
}

math::geom::Vector4D math::geom::operator*(const Vector4D& vec, Vector4D::value_type scalar) {
    return Vector4D(vec.x * scalar, vec.y * scalar, vec.z * scalar, vec.w * scalar);
}

math::geom::Vector4D math::geom::operator*(Vector4D::value_type scalar, const Vector4D& vec) {
    return vec * scalar;
}

math::geom::Vector4D math::geom::operator/(const Vector4D& vec, Vector4D::value_type scalar) {
    if (scalar == 0) throw std::invalid_argument("Division by zero");
    return vec * (1. / scalar);
}