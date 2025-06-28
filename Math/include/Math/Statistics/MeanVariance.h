#pragma once

#include "math_export.hpp"

#include <vector>

namespace math::stat {

    double MATH_EXPORT mean(const std::vector<double>& data);
    double MATH_EXPORT variance(const std::vector<double>& data);
    double MATH_EXPORT covariance(const std::vector<double>& x, const std::vector<double>& y);
    double MATH_EXPORT correlation(const std::vector<double>& x, const std::vector<double>& y);

} // namespace math::stat

double math::stat::mean(const std::vector<double>& data) {
    // ...
    return 0.0;
}

double math::stat::variance(const std::vector<double>& data) {
    // ...
    return 0.0;
}

double math::stat::covariance(const std::vector<double>& x, const std::vector<double>& y) {
    // ...
    return 0.0;
}

double math::stat::correlation(const std::vector<double>& x, const std::vector<double>& y) {
    // ...
    return 0.0;
}