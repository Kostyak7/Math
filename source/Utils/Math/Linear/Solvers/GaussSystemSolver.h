#pragma once 

#include "Interfaces/ILinearSystemSolver.h"

namespace math::linal {

    class GaussLinearSystemSolver final : public ILinearSystemSolver {
    public:
        GaussLinearSystemSolver(const Params& params = {});
        FVector solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& x0 = {}) override;
    };

} // namespace math::linal