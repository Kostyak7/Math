#pragma once 

#include <functional>

namespace math::calcls {

	double derivative(const std::function<double(double)>& f, double x, double h = 1e-5);

} // namespace math::calcls
