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

        FVector v_new = std::move(r);
        FVector w_new = apply_precondition(v_new);
        double beta_new2 = v_new.dot(w_new);
        if (beta_new2 < 0.0) {
            throw std::runtime_error("Preconditioner in MINRES is not positive definite");
        }

        FVector v_old(n, 0.0);
        FVector v(n, 0.0);
        double r_norm2 = r_norm * r_norm;
        FVector w(n, 0.0);

        double beta_new = std::sqrt(beta_new2);
        const double beta_one = beta_new;

        double c = 1.0; // косинус вращения Гивенса
        double c_old = 1.0;
        double s = 0.0; // синус вращения Гивенса
        double s_old = 0.0;
        FVector p_oold(n, 0.0);
        FVector p_old(n, 0.0);
        FVector p = p_old;
        double eta = 1.0;

        // порог сходимости
        const double threshold2 = m_iterative_solving_params.tolerance * m_iterative_solving_params.tolerance * rhs_norm * rhs_norm;

        size_t iter = 0; 
        for (; iter < max_iteration_count; ++iter) {
            // Предобусловленный алгоритм Ланцоша
            const double beta = beta_new;
            v_old = v;
            v_new /= beta_new;
            w_new /= beta_new;
            v = v_new;
            w = w_new;
            v_new = matrix * w - beta * v_old;

            const double alpha = v_new.dot(w);
            v_new -= alpha * v;
            w_new = apply_precondition(v_new);

            beta_new2 = v_new.dot(w_new);
            if (beta_new2 < 0.0) {
                throw std::runtime_error("Preconditioner in MINRES is not positive definite");
            }
            beta_new = std::sqrt(beta_new2);

            // Вращение Гивенса
            const double r2 = s * alpha + c * c_old * beta;
            const double r3 = s_old * beta;
            const double r1_hat = c * alpha - c_old * s * beta;
            const double r1 = std::hypot(r1_hat, beta_new);

            c_old = c;
            s_old = s;
            c = r1_hat / r1;
            s = beta_new / r1;

            p_oold = p_old;
            p_old = p;
            p = (w - r2 * p_old - r3 * p_oold) / r1;
            x += beta_one * c * eta * p;

            r_norm2 *= s * s;
            if (r_norm2 < threshold2) {
                break;
            }

            eta = -s * eta;
        }

        collect_solution_stats({ std::sqrt(r_norm2), iter});
        return x;
        }, matrix);
}
