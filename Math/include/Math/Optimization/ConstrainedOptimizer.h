#pragma once

#include "math_export.hpp"

#include <functional>
#include <vector>

namespace math::optim {

    class MATH_EXPORT ConstrainedOptimizer {
    public:
        static double optimize(const std::function<double(std::vector<double>)>& f,
            const std::vector<std::function<bool(std::vector<double>)>>& constraints,
            const std::vector<double>& init);
    };

} // namespace math::optim

double math::optim::ConstrainedOptimizer::optimize(
    const std::function<double(std::vector<double>)>& f,
    const std::vector<std::function<bool(std::vector<double>)>>& constraints,
    const std::vector<double>& init)
{
    // ...
    return 0.0;
}