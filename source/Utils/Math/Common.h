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

} // namespace math
