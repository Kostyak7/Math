#pragma once 

#include "Interfaces/IKrylovTypeLinearSystemSolver.h"

namespace math::linal {

    class CGLinearSystemSolver final : public IKrylovTypeLinearSystemSolver {
    public:
        CGLinearSystemSolver(size_t Krylov_subspace_dimension = 30, 
                             const IterativeSolvingParams& iter_params = {},
                             const Params& params = {}, 
                             std::unique_ptr<IPreconditioner> preconditioner = nullptr);
        DVector solve(const AnyMatrixConstRef& matrix, const DVector& rhs, const DVector& x0 = {}) override;

    private:
        bool slae_check(const IMatrix& matrix, const DVector& rhs) const override;
    };

} // namespace math::linal    
