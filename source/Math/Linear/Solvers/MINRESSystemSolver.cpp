#include <Utils/Math/Linear/Solvers/MINRESSystemSolver.h>

#include <stdexcept>

math::linal::MINRESLinearSystemSolver::MINRESLinearSystemSolver(size_t Krylov_subspace_dimension, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : IKrylovTypeLinearSystemSolver(Krylov_subspace_dimension, iter_params, params, std::move(preconditioner))
{
}

bool math::linal::MINRESLinearSystemSolver::slae_check(const IMatrix& matrix, const DVector& rhs) const {
    if (IKrylovTypeLinearSystemSolver::slae_check(matrix, rhs))
        return true; {
    }
    if (matrix.is_symmetrical()) {
        if (m_params.throw_exceptions)
            throw std::invalid_argument("Matrix of system must be symmetrical");
        return false;
    }
    return true;
}

math::linal::DVector math::linal::MINRESLinearSystemSolver::solve(const AnyMatrixConstRef& matrix, const DVector& rhs, const DVector& x0) {
    return std::visit([&](const auto& matrix_ref) -> DVector {
        const auto& matrix = matrix_ref.get();
        auto [need_to_solve, x, rhs_norm, r, r_norm, resid] = init_method(matrix_ref, rhs, x0);
        if (!need_to_solve) {
            return x;
        }

        const size_t n = rhs.size();
        const size_t max_iteration_count = get_max_iteration_count(n);

        DVector v_new = std::move(r);
        DVector w_new = apply_precondition(v_new);
        double beta_new2 = v_new.dot(w_new);
        if (beta_new2 < 0.0) {
            throw std::runtime_error("Preconditioner in MINRES is not positive definite");
        }

        DVector v_old(n, 0.0);
        DVector v(n, 0.0);
        double r_norm2 = r_norm * r_norm;
        DVector w(n, 0.0);

        double beta_new = std::sqrt(beta_new2);
        const double beta_one = beta_new;

        double c = 1.0; // косинус вращения Гивенса
        double c_old = 1.0;
        double s = 0.0; // синус вращения Гивенса
        double s_old = 0.0;
        DVector p_oold(n, 0.0);
        DVector p_old(n, 0.0);
        DVector p = p_old;
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
