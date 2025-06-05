#pragma once 

#include "IIterativeLinearSystemSolver.h"

namespace math::linal {

    class MATH_EXPORT IKrylovTypeLinearSystemSolver : public IIterativeLinearSystemSolver {
    public:
        IKrylovTypeLinearSystemSolver(size_t Krylov_subspace_dimension = 30, 
                                      const IterativeSolvingParams& iter_params = {}, 
                                      const Params& params = {}, 
                                      std::unique_ptr<IPreconditioner> preconditioner = nullptr);

    protected:
        size_t get_Krylov_subspace_dimension() const;
        size_t get_Krylov_subspace_dimension(size_t system_size) const;

    private:
        const size_t m_Krylov_subspace_dimension;
    };

} // namespace math::linal    
