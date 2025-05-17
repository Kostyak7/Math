#include "CGSystemSolver.h"

#include <stdexcept>

math::linal::CGLinearSystemSolver::CGLinearSystemSolver(size_t Krylov_subspace_dimension, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : IKrylovTypeLinearSystemSolver(Krylov_subspace_dimension, iter_params, params, std::move(preconditioner))
{
}

bool math::linal::CGLinearSystemSolver::slae_check(const IMatrix& matrix, const DVector& rhs) const {
    if (IKrylovTypeLinearSystemSolver::slae_check(matrix, rhs)) {
        return true;
    }
    if (matrix.is_symmetrical()) {
        if (m_params.throw_exceptions)
            throw std::invalid_argument("Matrix of system must be symmetrical");
        return false;
    }
    return true;
}

math::linal::DVector math::linal::CGLinearSystemSolver::solve(const AnyMatrixConstRef& matrix, const DVector& rhs, const DVector& x0) {
    return std::visit([&](const auto& matrix_ref) -> DVector {
        const auto& matrix = matrix_ref.get();
        auto [need_to_solve, x, rhs_norm, r, beta, resid] = init_method(matrix, rhs, x0);
        if (!need_to_solve) {
            return x;
        }

        const size_t n = rhs.size();
        const size_t max_iteration_count = get_max_iteration_count(n);

        DVector p(n);
        DVector q(n);
        DVector rhat(n);

        double rho, rho_prev = 0.0;

        for (size_t iter = 1; iter <= max_iteration_count; ++iter) {
            rhat = apply_precondition(r);

            rho = r.dot(rhat);

            if (iter == 1) {
                p = rhat;
            }
            else {
                p = rhat + p * (rho / rho_prev);
            }

            q = matrix * p;

            double alpha = rho / p.dot(q);

            x += p * alpha;
            r -= q * alpha;

            beta = r.norm();
            resid = beta / rhs_norm;
            if (check_convergence_criterion(r, rhs_norm, iter)) {
                collect_solution_stats({ resid, iter });
                return x;
            }

            rho_prev = rho;
        }

        collect_solution_stats({ resid, max_iteration_count });
        return x;
        }, matrix);
}
