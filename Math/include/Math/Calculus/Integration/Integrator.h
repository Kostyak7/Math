#pragma once 

#include "math_export.hpp"

#include <functional>

namespace math::calcls {

    class Integrator {
    public:
        static double MATH_EXPORT trapezoidal(const std::function<double(double)>& f, double a, double b, int n);
        static double MATH_EXPORT simpson(const std::function<double(double)>& f, double a, double b, int n);
    };

} // namespace math::calcls
