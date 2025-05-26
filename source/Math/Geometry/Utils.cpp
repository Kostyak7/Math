#include <Utils/Math/Geometry/Utils.h>

#include <stdexcept>

math::geom::Matrix2x2 math::geom::rotation_matrix_2d(double angle_rad) {
	Matrix2x2 result = identity_matrix<2>();
	double c = std::cos(angle_rad);
	double s = std::sin(angle_rad);

	result[0][0] = c;  result[0][1] = -s;
	result[1][0] = s;  result[1][1] = c;

	return result;
}

math::geom::Matrix3x3 math::geom::rotation_matrix_x(double angle_rad) {
	Matrix3x3 result = identity_matrix<3>();
	double c = std::cos(angle_rad);
	double s = std::sin(angle_rad);

	result[1][1] = c;  result[1][2] = -s;
	result[2][1] = s;  result[2][2] = c;

	return result;
}

math::geom::Matrix3x3 math::geom::rotation_matrix_y(double angle_rad) {
	Matrix3x3 result = identity_matrix<3>();
	double c = std::cos(angle_rad);
	double s = std::sin(angle_rad);

	result[1][1] = c;  result[1][2] = -s;
	result[2][1] = s;  result[2][2] = c;

	return result;
}

math::geom::Matrix3x3 math::geom::rotation_matrix_z(double angle_rad) {
	Matrix3x3 result = identity_matrix<3>();
	double c = std::cos(angle_rad);
	double s = std::sin(angle_rad);

	result[1][1] = c;  result[1][2] = -s;
	result[2][1] = s;  result[2][2] = c;

	return result;
}

math::geom::Vector2D math::geom::rotate_2d(const Vector2D& vec, double angle) {
	double cos_a = std::cos(angle);
	double sin_a = std::sin(angle);
	return {
		vec.x() * cos_a - vec.y() * sin_a,
		vec.x() * sin_a + vec.y() * cos_a
	};
}

math::geom::Vector3D math::geom::rotate_x(const Vector3D& vec, double angle) {
	double cos_a = std::cos(angle);
	double sin_a = std::sin(angle);
	return {
		vec.x(),
		vec.y() * cos_a - vec.z() * sin_a,
		vec.y() * sin_a + vec.z() * cos_a
	};
}

math::geom::Vector3D math::geom::rotate_y(const Vector3D& vec, double angle) {
	double cos_a = std::cos(angle);
	double sin_a = std::sin(angle);
	return {
		vec.x() * cos_a + vec.z() * sin_a,
		vec.y(),
		-vec.x() * sin_a + vec.z() * cos_a
	};
}

math::geom::Vector3D math::geom::rotate_z(const Vector3D& vec, double angle) {
	double cos_a = std::cos(angle);
	double sin_a = std::sin(angle);
	return {
		vec.x() * cos_a - vec.y() * sin_a,
		vec.x() * sin_a + vec.y() * cos_a,
		vec.z()
	};
}

math::geom::Point2D math::geom::rotate_around(const Point2D& p, const Point2D& center, double angle) {
	Vector2D rel = p - center;
	Vector2D rotated = rotate_2d(rel, angle);
	return center + rotated.get_end_point();
}

math::geom::Point3D math::geom::rotate_around_axis(const Point3D& p, const Vector3D& axis, double angle) {
	Vector3D u = normalized(axis);
	double cosA = std::cos(angle);
	double sinA = std::sin(angle);

	Vector3D v = p;
	return (u * (u.dot(v)) * (1 - cosA) + v * cosA + cross_product(u, v) * sinA).get_end_point();
}

math::geom::Point2D math::geom::project(const Point2D& point, const Line2D& line) {
	Vector2D ab(line.p1, line.p2);
	if (dcmp(ab.length()) == 0)
		return point;
	Vector2D ap(line.p1, point);
	double t = ap.dot(ab) / ab.dot(ab);
	return line.p1 + ab.get_end_point() * t;
}

math::geom::Matrix4x4 math::geom::orthographic_projection(double left, double right, double bottom, double top, double near, double far) {
	Matrix4x4 result = identity_matrix<4>();

	result[0][0] = 2.0f / (right - left);
	result[1][1] = 2.0f / (top - bottom);
	result[2][2] = -2.0f / (far - near);

	result[0][3] = -(right + left) / (right - left);
	result[1][3] = -(top + bottom) / (top - bottom);
	result[2][3] = -(far + near) / (far - near);

	return result;
}

math::geom::Matrix4x4 math::geom::perspective_projection(double fov_rad, double aspect, double near, double far) {
	Matrix4x4 result{};
	double tan_half_fov = std::tan(fov_rad / 2.0f);

	result[0][0] = 1.0f / (aspect * tan_half_fov);
	result[1][1] = 1.0f / tan_half_fov;
	result[2][2] = -(far + near) / (far - near);
	result[2][3] = -1.0f;
	result[3][2] = -(2.0f * far * near) / (far - near);

	return result;
}

math::geom::Matrix2x2 math::geom::scaling_matrix(double sx, double sy) {
	Matrix2x2 result{};
	result[0][0] = sx;
	result[1][1] = sy;
	return result;
}

math::geom::Matrix3x3 math::geom::scaling_matrix(double sx, double sy, double sz) {
	Matrix3x3 result{};
	result[0][0] = sx;
	result[1][1] = sy;
	result[2][2] = sz;
	return result;
}

math::geom::Matrix4x4 math::geom::scaling_matrix(double sx, double sy, double sz, double sw) {
	Matrix4x4 result{};
	result[0][0] = sx;
	result[1][1] = sy;
	result[2][2] = sz;
	result[3][3] = sw;
	return result;
}

math::geom::Matrix4x4 math::geom::scale_around_point(double sx, double sy, double sz, double cx, double cy, double cz) {
	Matrix4x4 translateToOrigin = translation_matrix(-cx, -cy, -cz);
	Matrix4x4 scale = scaling_matrix(sx, sy, sz, 1.0f);
	Matrix4x4 translateBack = translation_matrix(cx, cy, cz);
	return translateBack * scale * translateToOrigin;
}

math::geom::Matrix3x3 math::geom::translation_matrix(double x, double y) {
	Matrix3x3 result = identity_matrix<3>();
	result[0][2] = x;
	result[1][2] = y;
	return result;
}

math::geom::Matrix4x4 math::geom::translation_matrix(double x, double y, double z) {
	Matrix4x4 result = identity_matrix<4>();
	result[0][3] = x;
	result[1][3] = y;
	result[2][3] = z;
	return result;
}

math::geom::Vector2D math::geom::lerp(const Vector2D& a, const Vector2D& b, double t) {
	return a + (b - a) * t;
}

math::geom::Vector3D math::geom::lerp(const Vector3D& a, const Vector3D& b, double t) {
	return a + (b - a) * t;
}

double math::geom::signed_angle_2d(const Vector2D& v1, const Vector2D& v2) {
	double angle = angle_between(v1, v2);
	return (pseudodot(v1, v2) < 0) ? -angle : angle;
}

math::geom::Matrix4x4 math::geom::look_at(const Point3D& eye, const Point3D& target, const Vector3D& up) {
	Vector3D f = normalized(Vector3D{ target - eye });
	Vector3D r = normalized(cross_product(f, up));
	Vector3D u = cross_product(r, f);

	return {
		{	r.x(),	r.y(),	r.z(),	-r.dot(eye) },
		{	u.x(),	u.y(),	u.z(),	-u.dot(eye) },
		{  -f.x(), -f.y(), -f.z(),	 f.dot(eye) },
		{	  0.0,    0.0,    0.0,			1.0 },
	};
}

double math::geom::triangle_area(const Point2D& a, const Point2D& b, const Point2D& c) {
	return 0.5 * ((b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x()));
}

double math::geom::tetrahedron_volume(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4) {
	Vector3D v1 = p2 - p1;
	Vector3D v2 = p3 - p1;
	Vector3D v3 = p4 - p1;
	return v1.dot(cross_product(v2, v3)) / 6.0;
}

double math::geom::orientation_2d(const Point2D& a, const Point2D& b, const Point2D& c) {
	return (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x());
}

double math::geom::orientation_3d(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
	Vector3D ab = b - a;
	Vector3D ac = c - a;
	Vector3D ad = d - a;
	return ab.dot(cross_product(ac, ad));
}
