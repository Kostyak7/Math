#pragma once 

#include <Utils/Math/Linear/FVector.h>

namespace math::linal {

    class IConvergenceCriterion {
    public:
        virtual ~IConvergenceCriterion() = default;

        virtual bool is_converged(const FVector& residual, double rhs_norm, size_t iteration) const = 0;
    };

} // namespace math::linal