#pragma once 

#include <functional>

namespace math::calcls {

    class IODESolver {
    public:
        virtual double solve(const std::function<double(double, double)>& f, double y0, double x0, double x, double h = 1e-6) = 0;
    };

} // namespace math::calcls
