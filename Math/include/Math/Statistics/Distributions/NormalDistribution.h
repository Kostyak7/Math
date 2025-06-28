#pragma once

#include "math_export.hpp"

namespace math::stat {

    class MATH_EXPORT NormalDistribution {
    public:
        static double pdf(double x, double mean = 0.0, double stddev = 1.0);
        static double cdf(double x, double mean = 0.0, double stddev = 1.0);
    };

} // namespace math::stat

double math::stat::NormalDistribution::pdf(double x, double mean = 0.0, double stddev = 1.0) {
    // ...
    return 0.0;
}

double math::stat::NormalDistribution::cdf(double x, double mean = 0.0, double stddev = 1.0) {
    // ...
    return 0.0;
}