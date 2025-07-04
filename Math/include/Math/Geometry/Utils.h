#pragma once 

#include "Line.h"
#include "Matrix.h"
#include "Plane.h"

#include <Math/Common.h>

namespace math::geom {

    Matrix2x2 MATH_EXPORT rotation_matrix_2d(double angle_rad);
    Matrix3x3 MATH_EXPORT rotation_matrix_x(double angle_rad);
    Matrix3x3 MATH_EXPORT rotation_matrix_y(double angle_rad);
    Matrix3x3 MATH_EXPORT rotation_matrix_z(double angle_rad);

    Vector2D MATH_EXPORT rotate_2d(const Vector2D& vec, double angle);
    Vector3D MATH_EXPORT rotate_x(const Vector3D& vec, double angle);
    Vector3D MATH_EXPORT rotate_y(const Vector3D& vec, double angle);
    Vector3D MATH_EXPORT rotate_z(const Vector3D& vec, double angle);

    Point2D MATH_EXPORT rotate_around(const Point2D& p, const Point2D& center, double angle);
    Point3D MATH_EXPORT rotate_around_axis(const Point3D& p, const Vector3D& axis, double angle);

    Point2D MATH_EXPORT project(const Point2D& point, const Line2D& line);

    template <size_t N>
    Point<N> project_onto(const Point<N>& point, const Line<N>& line);

    Matrix4x4 MATH_EXPORT orthographic_projection(double left, double right, double bottom, double top, double near, double far);
    Matrix4x4 MATH_EXPORT perspective_projection(double fov_rad, double aspect, double near, double far);

    Matrix2x2 MATH_EXPORT scaling_matrix(double sx, double sy);
    Matrix3x3 MATH_EXPORT scaling_matrix(double sx, double sy, double sz);
    Matrix4x4 MATH_EXPORT scaling_matrix(double sx, double sy, double sz, double sw);

    Matrix4x4 MATH_EXPORT scale_around_point(double sx, double sy, double sz, double cx, double cy, double cz);

    Matrix3x3 MATH_EXPORT translation_matrix(double x, double y);
    Matrix4x4 MATH_EXPORT translation_matrix(double x, double y, double z);

    Vector2D MATH_EXPORT lerp(const Vector2D& a, const Vector2D& b, double t);
    Vector3D MATH_EXPORT lerp(const Vector3D& a, const Vector3D& b, double t);

    template <size_t N>
    double angle_between(const Vector<N>& v1, const Vector<N>& v2);

    double MATH_EXPORT signed_angle_2d(const Vector2D& v1, const Vector2D& v2);

    template <size_t N>
    double distance_to_line_nd(const Point<N>& point, const Line<N>& line);
    template <size_t N>
    bool is_point_on_line(const Point<N>& p, const Line<N>& line, double eps = TOLERANCE);

    Matrix4x4 MATH_EXPORT look_at(const Point3D& eye, const Point3D& target, const Vector3D& up);

    double MATH_EXPORT triangle_area(const Point2D& a, const Point2D& b, const Point2D& c);
    double MATH_EXPORT tetrahedron_volume(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4);

    template <size_t N>
    Point<N> closest_point_on_line(const Point<N>& point, const Line<N>& line);

    double MATH_EXPORT orientation_2d(const Point2D& a, const Point2D& b, const Point2D& c); // > 0 ����� �������, < 0 ������ �������, = 0 �� ����� ������
    double MATH_EXPORT orientation_3d(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d); // > 0 �� ���� �� ������ ��������� ABC, ..., = 0 ��� � ����� ���������

} // namespace math::geom

#include "Utils.tpp"
