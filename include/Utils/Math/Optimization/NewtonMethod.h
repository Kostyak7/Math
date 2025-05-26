#pragma once

#include <functional>

namespace math::optim {

    class NewtonMethod {
    public:
        static double optimize(const std::function<double(double)>& f, const std::function<double(double)>& df, const std::function<double(double)>& ddf, double init, int steps);
    };

} // namespace math::optim

double math::optim::NewtonMethod::optimize(const std::function<double(double)>& f, const std::function<double(double)>& df, const std::function<double(double)>& ddf, double init, int steps) {
    // ...
    return 0.0;
}