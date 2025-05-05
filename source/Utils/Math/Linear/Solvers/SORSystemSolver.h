#pragma once 

#include "Interfaces/IIterativeLinearSystemSolver.h"

namespace math::linal {

    class SORLinearSystemSolver final : public IIterativeLinearSystemSolver {
    public:
        SORLinearSystemSolver(double w = 1.0, 
                              size_t auto_tune_period = 10, 
                              const IterativeSolvingParams& iter_params = {},
                              const Params& params = {},
                              std::unique_ptr<IPreconditioner> preconditioner = nullptr);
        FVector solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& x0 = {}) override;

    private:
        void adapt_relaxation_parameter(double residual);

    private:
        double m_w;
        size_t m_auto_tune_period;
    };

} // namespace math::linal    
