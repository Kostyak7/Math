#pragma once 

#include "GMRESSystemSolver.h"

namespace math::linal {

    class FGMRESLinearSystemSolver final : public GMRESLinearSystemSolver {
    public:
        FGMRESLinearSystemSolver(size_t Krylov_subspace_dimension = 30, 
                                 const IterativeSolvingParams& iter_params = {}, 
                                 const Params& params = {}, 
                                 std::unique_ptr<IPreconditioner> preconditioner = nullptr);
        FVector solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& x0 = {}) override;

    private:
        void compute_correction(size_t k, size_t n, const std::vector<double>& H, size_t ldH, const std::vector<double>& s, const std::vector<FVector>& Z, FVector& x) const override;
    };

} // namespace math::linal    
