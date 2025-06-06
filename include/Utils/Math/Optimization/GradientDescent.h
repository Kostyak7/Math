#pragma once

#include "math_export.hpp"

#include <functional>

namespace math::optim {

    class MATH_EXPORT GradientDescent {
    public:
        static double optimize(const std::function<double(double)>& f, const std::function<double(double)>& df, double init, double lr, int steps);
    };

} // namespace math::optim

double math::optim::GradientDescent::optimize(const std::function<double(double)>& f, const std::function<double(double)>& df, double init, double lr, int steps) {
    // ...
    return 0.0;
}
