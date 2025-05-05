#include "BiCGSTABSystemSolver.h"

#include <stdexcept>

math::linal::BiCGSTABLinearSystemSolver::BiCGSTABLinearSystemSolver(size_t Krylov_subspace_dimension, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : IKrylovTypeLinearSystemSolver(Krylov_subspace_dimension, iter_params, params, std::move(preconditioner))
{
}
    
math::linal::FVector math::linal::BiCGSTABLinearSystemSolver::solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& x0) {
    return std::visit([this, &rhs, &x0](const auto& matrix) -> FVector {
        auto [need_to_solve, x, rhs_norm, r, r_norm, resid] = init_method(matrix, rhs, x0);
        if (!need_to_solve) {
            return x;
        }

        const size_t n = rhs.size();
        const size_t max_iteration_count = get_max_iteration_count(n);

        FVector r0 = r;
        FVector v(n);
        FVector p(n);
        FVector phat(n);
        FVector s(n);
        FVector shat(n);
        FVector t(n);

        double alpha = 0.0;
        double omega = 0.0;
        double rho_prev = 0.0;

        for (size_t iter = 1; iter <= max_iteration_count; ++iter) {
            double rho = r0.dot(r);
            if (dcmp(rho) == 0) {
                collect_solution_stats({ r.norm() / rhs_norm, iter });
                return x; // ������: ��������� ������������ ������� ����
            }

            // ���������� ����������� p
            if (iter == 1) {
                p = r; 
            }
            else {
                double beta = (rho / rho_prev) * (alpha / omega);
                p = r + (p - v * omega) * beta;
            }

            // ���������� �������������������: phat = C*p
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

            // ���������� ������������������� � s: shat = C*s
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
                return x; // ������: �������� omega ������� ���
            }

            rho_prev = rho;
        }

        collect_solution_stats({ resid, max_iteration_count });
        return x; // ��������� ������������ ����� ��������
        }, matrix);
}
