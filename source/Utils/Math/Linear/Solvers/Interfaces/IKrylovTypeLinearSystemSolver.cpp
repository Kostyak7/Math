#include "IKrylovTypeLinearSystemSolver.h"

math::linal::IKrylovTypeLinearSystemSolver::IKrylovTypeLinearSystemSolver(size_t Krylov_subspace_dimension, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : IIterativeLinearSystemSolver(iter_params, params, std::move(preconditioner))
    , m_Krylov_subspace_dimension(Krylov_subspace_dimension)
{
}

size_t math::linal::IKrylovTypeLinearSystemSolver::get_Krylov_subspace_dimension() const {
    return m_Krylov_subspace_dimension;
}

size_t math::linal::IKrylovTypeLinearSystemSolver::get_Krylov_subspace_dimension(size_t system_size) const {
    return std::min(m_Krylov_subspace_dimension, system_size);
}