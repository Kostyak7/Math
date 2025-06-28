#pragma once 

#include "math_export.hpp"

#include <functional>

namespace math::calcls {

	double MATH_EXPORT derivative(const std::function<double(double)>& f, double x, double h = 1e-5);

} // namespace math::calcls
