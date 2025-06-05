#pragma once

#include "math_export.hpp"

#include <vector>

namespace math::stat {

    class MATH_EXPORT Histogram {
    public:
        static std::vector<int> build(const std::vector<double>& data, int bins);
    };

} // namespace math::stat

std::vector<int> math::stat::Histogram::build(const std::vector<double>&data, int bins) {
    // ...
    return {};
}

