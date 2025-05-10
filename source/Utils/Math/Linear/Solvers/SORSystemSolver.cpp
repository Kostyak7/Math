#include "SORSystemSolver.h"

#include <stdexcept>

namespace {

    void do_sor_step(const math::linal::BandMatrix& matrix, math::linal::DVector& x, const math::linal::DVector& rhs, double w) {
        const size_t n = matrix.get_height();
        const size_t isl = matrix.get_isl();

        double sigma;

        for (size_t i = 0; i < n; ++i) {
            sigma = 0.0;
            for (size_t j = (i < isl) ? 0 : i - isl + 1; j < i; ++j) {
                sigma += matrix[j][i - j] * x[j];
            }

            const auto& row = matrix[i];
            for (size_t k = 1; k < isl; ++k) {
                size_t j = i + k;
                if (j < n) {
                    sigma += row[k] * x[j];
                }
            }
            x[i] = (1 - w) * x[i] + (w / matrix[i][0]) * (rhs[i] - sigma);
        }
    }

    void do_sor_step(const math::linal::DenseMatrix& matrix, math::linal::DVector& x, const math::linal::DVector& rhs, double w) {
        const size_t n = matrix.get_height();

        for (size_t i = 0; i < n; ++i) {
            double sigma = 0.0;
            for (size_t j = 0; j < n; ++j) {
                if (j != i) {
                    sigma += matrix[i][j] * x[j];
                }
            }
            x[i] = (1 - w) * x[i] + (w / matrix[i][i]) * (rhs[i] - sigma);
        }
    }

    void do_sor_step(const math::linal::SparseMatrix& matrix, math::linal::DVector& x, const math::linal::DVector& rhs, double w) {
        const size_t n = matrix.get_height();

        // ...
    }

} // namespace 

math::linal::SORLinearSystemSolver::SORLinearSystemSolver(double w, size_t auto_tune_period, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : IIterativeLinearSystemSolver(iter_params, params, std::move(preconditioner))
    , m_w(w)
    , m_auto_tune_period(auto_tune_period)
{
    if (w <= 0 || w >= 2) {
        if (m_params.throw_exceptions)
            throw std::invalid_argument("Relaxation parameter w must be in (0, 2)");
        else
            m_w = 1.0;
    }
}

math::linal::DVector math::linal::SORLinearSystemSolver::solve(const AnyMatrix& matrix, const DVector& rhs, const DVector& x0) {
    return std::visit([this, &rhs, &x0](const auto& matrix) -> DVector {
        auto [need_to_solve, x, rhs_norm, r, r_norm, resid] = init_method(matrix, rhs, x0);
        if (!need_to_solve) {
            return x;
        }

        const size_t n = rhs.size();
        const size_t max_iteration_count = get_max_iteration_count(n);

        DVector z;

        size_t iter = 0;
        do {
            z = apply_precondition(r);

            do_sor_step(matrix, x, z, m_w);

            // Автоматическая адаптация параметра релаксации
            if (m_auto_tune_period != 0 && iter % m_auto_tune_period == 0) {
                adapt_relaxation_parameter(resid);
            }
            ++iter;

            r = rhs - matrix * x;
            resid = r.norm() / rhs_norm;

        } while (!check_convergence_criterion(r, rhs_norm, iter) && iter <= m_iterative_solving_params.max_iteration_count);

        collect_solution_stats({ m_iterative_solving_params.tolerance, iter });
        return x;
        }, matrix);
}

void math::linal::SORLinearSystemSolver::adapt_relaxation_parameter(double residual) {
    static double last_residual = std::numeric_limits<double>::max();

    if (residual >= last_residual * 0.99) { // чисто эвристика
        m_w = std::max(0.95, m_w * 0.95);
    }
    else {
        m_w = std::min(1.95, m_w * 1.01);
    }

    last_residual = residual;
}