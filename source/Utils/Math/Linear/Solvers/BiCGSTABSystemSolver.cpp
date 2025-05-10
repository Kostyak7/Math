#include "BiCGSTABSystemSolver.h"

#include <stdexcept>

math::linal::BiCGSTABLinearSystemSolver::BiCGSTABLinearSystemSolver(size_t Krylov_subspace_dimension, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : IKrylovTypeLinearSystemSolver(Krylov_subspace_dimension, iter_params, params, std::move(preconditioner))
{
}
    
math::linal::DVector math::linal::BiCGSTABLinearSystemSolver::solve(const AnyMatrix& matrix, const DVector& rhs, const DVector& x0) {
    return std::visit([this, &rhs, &x0](const auto& matrix) -> DVector {
        auto [need_to_solve, x, rhs_norm, r, r_norm, resid] = init_method(matrix, rhs, x0);
        if (!need_to_solve) {
            return x;
        }

        const size_t n = rhs.size();
        const size_t max_iteration_count = get_max_iteration_count(n);

        DVector r0 = r;
        DVector v(n);
        DVector p(n);
        DVector phat(n);
        DVector s(n);
        DVector shat(n);
        DVector t(n);

        double alpha = 0.0;
        double omega = 0.0;
        double rho_prev = 0.0;

        for (size_t iter = 1; iter <= max_iteration_count; ++iter) {
            double rho = r0.dot(r);
            if (dcmp(rho) == 0) {
                collect_solution_stats({ r.norm() / rhs_norm, iter });
                return x; // Ошибка: скалярное произведение слишком мало
            }

            // Вычисление направления p
            if (iter == 1) {
                p = r; 
            }
            else {
                double beta = (rho / rho_prev) * (alpha / omega);
                p = r + (p - v * omega) * beta;
            }

            // Применение предобуславливателя: phat = C*p
            phat = apply_precondition(p);

            v = matrix * phat;

            alpha = rho / r0.dot(v);

            s = r - v * alpha;

            resid = s.norm() / rhs_norm;
            if (check_convergence_criterion(s, rhs_norm, iter)) {
                x += phat * alpha;
                collect_solution_stats({ resid, iter });
                return x;
            }

            // Применение предобуславливателя к s: shat = C*s
            shat = apply_precondition(s);

            t = matrix * shat;

            omega = t.dot(s) / t.dot(t);

            x += phat * alpha + shat * omega;

            r = s - t * omega;

            resid = r.norm() / rhs_norm;
            if (check_convergence_criterion(r, rhs_norm, iter)) {
                collect_solution_stats({ resid, iter });
                return x;
            }

            if (dcmp(omega) == 0) {
                collect_solution_stats({ resid, iter });
                return x; // Ошибка: параметр omega слишком мал
            }

            rho_prev = rho;
        }

        collect_solution_stats({ resid, max_iteration_count });
        return x; // Превышено максимальное число итераций
        }, matrix);
}
