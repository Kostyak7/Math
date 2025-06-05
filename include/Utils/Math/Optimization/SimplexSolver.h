#pragma once

#include "math_export.hpp"

#include <vector>

namespace math::optim {

    class MATH_EXPORT SimplexSolver {
    public:
        static std::vector<double> solve(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const std::vector<double>& c);
    };

} // namespace math::optim

std::vector<double> math::optim::SimplexSolver::solve(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const std::vector<double>& c) {
    // ...
    return {};
}
