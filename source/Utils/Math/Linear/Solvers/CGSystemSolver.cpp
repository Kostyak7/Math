#include "CGSystemSolver.h"

#include <stdexcept>

math::linal::CGLinearSystemSolver::CGLinearSystemSolver(size_t Krylov_subspace_dimension, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : IKrylovTypeLinearSystemSolver(Krylov_subspace_dimension, iter_params, params, std::move(preconditioner))
{
}

math::linal::FVector math::linal::CGLinearSystemSolver::solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& x0) {
    return std::visit([this, &rhs, &x0](const auto& matrix) -> FVector {
        auto [need_to_solve, x, rhs_norm, r, beta, resid] = init_method(matrix, rhs, x0);
        if (!need_to_solve) {
            return x;
        }

        const size_t n = rhs.size();
        const size_t max_iteration_count = get_max_iteration_count(n);

        FVector p(n);
        FVector q(n);
        FVector rhat(n);

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
