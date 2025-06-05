#pragma once

#include "Constants.h"

#include <cmath>
#include <complex>

namespace math {

	constexpr auto TOLERANCE = 10e-10;

    inline int dcmp(double v1, double v2 = 0., double tolerance = TOLERANCE) {
        if (auto d = v1 - v2; fabs(d) > tolerance) {
            return (d > 0.) ? 1 : -1;
        }
        return 0;
    }

    template <typename T>
    int sign(T x) {
        return (x > T(0)) - (x < T(0));
    }

    inline double to_radians(double degrees) {
        return degrees * pi / 180.0;
    }

    inline double to_degrees(double radians) {
        return radians * 180.0 / pi;
    }

    inline double clamp(double value, double min_val, double max_val) {
        return std::max(min_val, std::min(max_val, value));
    }

    inline bool is_near(double a, double b, double eps = TOLERANCE) {
        return fabs(a - b) < eps;
    }

    inline double lerp(double a, double b, double t) {
        return a + (b - a) * t;
    }

} // namespace math
