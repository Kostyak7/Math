#pragma once 

#include "Interfaces/IKrylovTypeLinearSystemSolver.h"

namespace math::linal {

    class GMRESLinearSystemSolver : public IKrylovTypeLinearSystemSolver {
    public:
        GMRESLinearSystemSolver(size_t Krylov_subspace_dimension = 30, 
                                const IterativeSolvingParams& iter_params = {},
                                const Params& params = {}, 
                                std::unique_ptr<IPreconditioner> preconditioner = nullptr);
        DVector solve(const AnyMatrix& matrix, const DVector& rhs, const DVector& x0 = {}) override;

    protected:
        virtual void compute_correction(size_t k, size_t n, const std::vector<double>& H, size_t ldH, const std::vector<double>& s, const std::vector<DVector>& V, DVector& x) const;

        static DVector solve_H_system(size_t k, const std::vector<double>& H, size_t ldH, const std::vector<double>& s, bool throw_exception = true);
        static std::pair<double, double> generate_Givens_rotation(double dx, double dy);
        static void apply_Givens_rotation(double& dx, double& dy, double cs, double sn);
    };

} // namespace math::linal    
