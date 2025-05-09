#pragma once 

#include "Interfaces/IKrylovTypeLinearSystemSolver.h"

namespace math::linal {

    class CGLinearSystemSolver final : public IKrylovTypeLinearSystemSolver {
    public:
        CGLinearSystemSolver(size_t Krylov_subspace_dimension = 30, 
                             const IterativeSolvingParams& iter_params = {},
                             const Params& params = {}, 
                             std::unique_ptr<IPreconditioner> preconditioner = nullptr);
        FVector solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& x0 = {}) override;

    private:
        bool slae_check(const AnyMatrix& matrix, const FVector& rhs) const override;
    };

} // namespace math::linal    
