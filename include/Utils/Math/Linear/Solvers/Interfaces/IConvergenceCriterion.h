#pragma once 

#include <Utils/Math/Linear/DVector.h>

namespace math::linal {

    class MATH_EXPORT IConvergenceCriterion {
    public:
        virtual ~IConvergenceCriterion() = default;
        virtual bool is_converged(const DVector& residual, double rhs_norm, size_t iteration) const = 0;
    };

} // namespace math::linal