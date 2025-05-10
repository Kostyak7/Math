#pragma once 

#include "Interfaces/IKrylovTypeLinearSystemSolver.h"

namespace math::linal {

    class BiCGSTABLinearSystemSolver final : public IKrylovTypeLinearSystemSolver {
    public:
        BiCGSTABLinearSystemSolver(size_t Krylov_subspace_dimension = 30, 
                                   const IterativeSolvingParams& iter_params = {}, 
                                   const Params& params = {}, 
                                   std::unique_ptr<IPreconditioner> preconditioner = nullptr);
        DVector solve(const AnyMatrix& matrix, const DVector& rhs, const DVector& x0 = {}) override;
    };

} // namespace math::linal
