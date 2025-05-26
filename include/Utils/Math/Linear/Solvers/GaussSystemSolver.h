#pragma once 

#include "Interfaces/ILinearSystemSolver.h"

namespace math::linal {

    class GaussLinearSystemSolver final : public ILinearSystemSolver {
    public:
        GaussLinearSystemSolver(const Params& params = {});
        DVector solve(const AnyMatrixConstRef& matrix, const DVector& rhs, const DVector& x0 = {}) override;
    };

} // namespace math::linal