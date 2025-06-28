#pragma once 

#include "IODESolver.h"

namespace math::calcls {

    class MATH_EXPORT EulerSolver : public IODESolver {
    public:
        double solve(const std::function<double(double, double)>& f, double y0, double x0, double x, double h = 1e-6) override;
    };

} // namespace math::calcls
