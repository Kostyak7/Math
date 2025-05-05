#pragma once 

#include "Point.h"

#include <vector>

namespace math::geom {

	Vector2D rotate(const Vector2D& vec, double angle);

	Point2D project(const Point2D& point, const Line& line);

}