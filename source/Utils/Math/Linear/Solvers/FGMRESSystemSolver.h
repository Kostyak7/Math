#pragma once 

#include "GMRESSystemSolver.h"

namespace math::linal {

    class FGMRESLinearSystemSolver final : public GMRESLinearSystemSolver {
    public:
        FGMRESLinearSystemSolver(size_t Krylov_subspace_dimension = 30, 
                                 const IterativeSolvingParams& iter_params = {}, 
                                 const Params& params = {}, 
                                 std::unique_ptr<IPreconditioner> preconditioner = nullptr);
        DVector solve(const AnyMatrix& matrix, const DVector& rhs, const DVector& x0 = {}) override;

    private:
        void compute_correction(size_t k, size_t n, const std::vector<double>& H, size_t ldH, const std::vector<double>& s, const std::vector<DVector>& Z, DVector& x) const override;
    };

} // namespace math::linal    
