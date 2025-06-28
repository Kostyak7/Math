#pragma once

#include "math_export.hpp"

#include <vector>

namespace math::stat {

    class MATH_EXPORT RandomSampling {
    public:
        static std::vector<double> uniform(int n, double a, double b);
        static std::vector<double> normal(int n, double mean = 0.0, double stddev = 1.0);
    };

} // namespace math::stat

std::vector<double> math::stat::RandomSampling::uniform(int n, double a, double b) {
    // ...
    return {};
}

std::vector<double> math::stat::RandomSampling::normal(int n, double mean = 0.0, double stddev = 1.0) {
    // ...
    return {};
}