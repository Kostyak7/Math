#pragma once 

#include "Point.h"

#include <vector>

namespace math::geom {

	std::vector<Point2D> Graham_scan(std::vector<Point2D> points);

	std::vector<Point2D> Andrew_scan(std::vector<Point2D> points);

}