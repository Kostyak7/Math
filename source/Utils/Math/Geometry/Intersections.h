#pragma once 

#include "Line.h"
#include "Plane.h"

#include <optional>

namespace math::geom {

	std::optional<Point2D> segment_intersection(const Line& l1, const Line& l2);

}
