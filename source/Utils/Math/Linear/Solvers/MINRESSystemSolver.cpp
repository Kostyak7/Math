#include "MINRESSystemSolver.h"

#include <stdexcept>

math::linal::MINRESLinearSystemSolver::MINRESLinearSystemSolver(size_t Krylov_subspace_dimension, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : IKrylovTypeLinearSystemSolver(Krylov_subspace_dimension, iter_params, params, std::move(preconditioner))
{
}

math::linal::FVector math::linal::MINRESLinearSystemSolver::solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& x0) {
    return std::visit([this, &rhs, &x0](const auto& matrix) -> FVector {
        auto [need_to_solve, x, rhs_norm, r, r_norm, resid] = init_method(matrix, rhs, x0);
        if (!need_to_solve) {
            return x;
        }

        const size_t n = rhs.size();
        const size_t max_iteration_count = get_max_iteration_count(n);

        auto z = apply_precondition(r);
        double beta = sqrt(r.dot(z));
        if (beta <  m_iterative_solving_params.tolerance) {
            collect_solution_stats({ resid, 0 });
            return x;
        }

        FVector p0 = r;
        FVector s0 = matrix * p0;
        FVector p1 = p0;
        FVector s1 = s0;
        FVector p2, s2;

        size_t iter = 0;
        for (; iter < max_iteration_count; ++iter) {
            p2 = p1; p1 = p0;
            s2 = s1; s1 = s0;

            double alpha = r.dot(s1) / s1.dot(s1);
            x += alpha * p1;
            r -= alpha * s1;

            resid = r.norm() / rhs_norm;
            if (check_convergence_criterion(r, rhs_norm, iter)) {
                break;
            }

            p0 = s1;
            s0 = matrix * s1;

            double beta1 = s0.dot(s1) / s1.dot(s1);
            p0 = p0 - beta1 * p1;
            s0 = s0 - beta1 * s1;

            if (iter > 0) {
                double beta2 = s0.dot(s2) / s2.dot(s2);
                p0 = p0 - beta2 * p2;
                s0 = s0 - beta2 * s2;
            }
        }

        collect_solution_stats({ resid, iter });
        return x;
        }, matrix);
}
