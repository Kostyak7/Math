#include <Utils/Math/Linear/Solvers/Interfaces/IIterativeLinearSystemSolver.h>

math::linal::IIterativeLinearSystemSolver::IIterativeLinearSystemSolver(const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : ILinearSystemSolver(params)
    , m_iterative_solving_params(iter_params)
    , m_preconditioner(std::move(preconditioner))
{
}

math::linal::IIterativeLinearSystemSolver::IterativeSolvingParams math::linal::IIterativeLinearSystemSolver::get_iterative_solving_params() const {
    return m_iterative_solving_params;
}

math::linal::IIterativeLinearSystemSolver::SolutionStats math::linal::IIterativeLinearSystemSolver::get_last_solution_stats() const {
    return m_solution_stats;
}

math::linal::IIterativeLinearSystemSolver::Data math::linal::IIterativeLinearSystemSolver::init_method(const AnyMatrixConstRef& matrix, const DVector& rhs, const DVector& x0) {
    return std::visit([&](const auto& matrix_ref) -> Data {
        const auto& matrix = matrix_ref.get();
        Data res{ false, {}, 0.0, {}, 0.0, 0.0 };
        if (!slae_check(matrix, rhs)) {
            return res;
        }

        if (m_iterative_solving_params.additional_prechecks && matrix.is_identity()) {
            collect_solution_stats({ 0.0, 0 });
            res.x = rhs;
            return res;
        }

        const size_t n = rhs.size();

        res.rhs_norm = rhs.norm();
        if (m_iterative_solving_params.additional_prechecks && res.rhs_norm < m_iterative_solving_params.tolerance) {
            collect_solution_stats({ 0.0, 0 });
            res.x = DVector(n, 0.0);
            return res;
        }
        res.x = (x0.size() != n) ? DVector(n, 0.0) : x0;

        res.r = rhs - matrix * res.x;
        res.r_norm = res.r.norm();
        res.resid = res.r_norm / res.rhs_norm;
        if (check_convergence_criterion(res.r, res.rhs_norm, 0)) {
            collect_solution_stats({ res.resid, 0 });
            return res;
        }

        init_preconditioner(matrix_ref);

        res.need_to_solve = true;
        return res;
        }, matrix);
}

void math::linal::IIterativeLinearSystemSolver::init_preconditioner(const AnyMatrixConstRef& matrix) const {
    if (m_preconditioner) {
        m_preconditioner->init(matrix);
    }
}

math::linal::DVector math::linal::IIterativeLinearSystemSolver::apply_precondition(const DVector& vector) const {
    return m_preconditioner ? m_preconditioner->apply(vector) : vector;
}

bool math::linal::IIterativeLinearSystemSolver::check_convergence_criterion(const DVector& residual, double rhs_norm, size_t iteration) const {
    static auto default_convergence_criterion = std::make_unique<RelativeResidualCriterion>(m_iterative_solving_params.tolerance);
    if (m_iterative_solving_params.convergence_criterion) {
        return m_iterative_solving_params.convergence_criterion->is_converged(residual, rhs_norm, iteration);
    }
    return default_convergence_criterion->is_converged(residual, rhs_norm, iteration);
}

size_t math::linal::IIterativeLinearSystemSolver::get_max_iteration_count(size_t system_size) const {
    if (m_iterative_solving_params.max_iteration_count == 0) {
        return m_iterative_solving_params.coef_iteration_count * system_size;
    }
    return std::min(m_iterative_solving_params.max_iteration_count, m_iterative_solving_params.coef_iteration_count * system_size);
}

void math::linal::IIterativeLinearSystemSolver::collect_solution_stats(const SolutionStats& stats) {
    m_solution_stats = stats;
}