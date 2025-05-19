#include "FGMRESSystemSolver.h"

#include <stdexcept>

math::linal::FGMRESLinearSystemSolver::FGMRESLinearSystemSolver(size_t Krylov_subspace_dimension, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : GMRESLinearSystemSolver(Krylov_subspace_dimension, iter_params, params, std::move(preconditioner))
{
}

math::linal::DVector math::linal::FGMRESLinearSystemSolver::solve(const AnyMatrixConstRef& matrix, const DVector& rhs, const DVector& x0) {
    return std::visit([&](const auto& matrix_ref) -> DVector {
        const auto& matrix = matrix_ref.get();
        auto [need_to_solve, x, rhs_norm, r, beta, resid] = init_method(matrix_ref, rhs, x0);
        if (!need_to_solve) {
            return x;
        }

        const size_t n = rhs.size();
        const size_t m = get_Krylov_subspace_dimension(n);
        const size_t max_iteration_count = get_max_iteration_count(n);        

        std::vector<DVector> V(m + 1, DVector(n, 0.0));
        std::vector<DVector> Z(m + 1, DVector(n, 0.0));
        std::vector<double> H((m + 1) * m);
        std::vector<double> cs(m + 1);
        std::vector<double> sn(m + 1);
        std::vector<double> s(m + 1);

        for (size_t j = 1; j <= max_iteration_count;) {
            V[0] = r / beta; // первый ортонормированный вектор

            s[0] = beta;
            for (size_t i = 1; i <= m; ++i) {
                s[i] = 0.0;
            }

            for (size_t i = 0; i < m && j <= max_iteration_count; ++i, ++j) {
                Z[i] = apply_precondition(V[i]);

                DVector w = matrix * Z[i];

                // Ортогонализация w относительно предыдущих векторов V
                for (size_t k = 0; k <= i; ++k) {
                    H[k + i * (m + 1)] = w.dot(V[k]);
                    w -= V[k] * H[k + i * (m + 1)];
                }

                H[(i + 1) + i * (m + 1)] = w.norm();
                V[i + 1] = w * (1.0 / H[(i + 1) + i * (m + 1)]);

                // Применяем предыдущие вращения Гивенса к последнему столбцу H
                for (size_t k = 0; k < i; ++k) {
                    apply_Givens_rotation(H[k + i * (m + 1)], H[k + 1 + i * (m + 1)], cs[k], sn[k]);
                }

                // Генерируем новое вращение Гивенса
                std::tie(cs[i], sn[i]) = generate_Givens_rotation(H[i + i * (m + 1)], H[(i + 1) + i * (m + 1)]);

                // Применяем его к H и s
                apply_Givens_rotation(H[i + i * (m + 1)], H[(i + 1) + i * (m + 1)], cs[i], sn[i]);
                apply_Givens_rotation(s[i], s[i + 1], cs[i], sn[i]);

                resid = fabs(s[i + 1] / rhs_norm);
                if (resid < m_iterative_solving_params.tolerance) {
                    compute_correction(i + 1, n, H, m + 1, s, Z, x);
                    collect_solution_stats({ resid, j });
                    return x;
                }
            }

            if (m > 0) {
                compute_correction(m, n, H, m + 1, s, Z, x);
            }

            r = rhs - matrix * x;
            beta = r.norm();
            resid = beta / rhs_norm;
            if (check_convergence_criterion(r, rhs_norm, j)) {
                collect_solution_stats({ resid, j });
                return x;
            }
        }

        collect_solution_stats({ resid, max_iteration_count });
        return x;
        }, matrix);
}

// Решение верхней треугольной системы H y = s и обновление x += Z y
void math::linal::FGMRESLinearSystemSolver::compute_correction(size_t k, size_t n, const std::vector<double>& H, size_t ldH, const std::vector<double>& s, const std::vector<DVector>& Z, DVector& x) const {
    auto y = solve_H_system(k, H, ldH, s, m_params.throw_exceptions);

    for (size_t i = 0; i < k; ++i) {
        x += Z[i] * y[i];
    }
}
