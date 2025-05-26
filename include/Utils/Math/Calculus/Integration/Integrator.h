#pragma once 

#include <functional>

namespace math::calcls {

    class Integrator {
    public:
        static double trapezoidal(const std::function<double(double)>& f, double a, double b, int n);
        static double simpson(const std::function<double(double)>& f, double a, double b, int n);
    };

} // namespace math::calcls
