#pragma once 

#include "Point.h"

#include "math_export.hpp"

#include <vector>

namespace math::geom {

	std::vector<Point2D> MATH_EXPORT Graham_scan(std::vector<Point2D> points);

	std::vector<Point2D> MATH_EXPORT Andrew_scan(std::vector<Point2D> points);

}