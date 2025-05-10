#pragma once 

#include "ILinearSystemSolver.h"

#include <Utils/Math/Linear/Solvers/Preconditioners/IPreconditioner.h>
#include <Utils/Math/Linear/Solvers/ConvergenceCriterion.h>
#include <Utils/Math/Common.h>

#include <memory>

namespace math::linal {

    class IIterativeLinearSystemSolver : public ILinearSystemSolver {
    public:

        struct IterativeSolvingParams {
            double tolerance = TOLERANCE;
            size_t max_iteration_count = 40000;
            size_t coef_iteration_count = 7;
            bool additional_prechecks = true;
            std::shared_ptr<IConvergenceCriterion> convergence_criterion = std::make_shared<RelativeResidualCriterion>();
        };

        struct SolutionStats {
            double accuracy = 0.0;
            size_t iteration_count = 0;
        };

    public:
        IIterativeLinearSystemSolver(const IterativeSolvingParams& iter_params = {}, 
                                     const Params& params = {},
                                     std::unique_ptr<IPreconditioner> preconditioner = nullptr);

        IterativeSolvingParams get_iterative_solving_params() const;
        SolutionStats get_last_solution_stats() const;

    protected:
        struct Data {
            bool need_to_solve;
            DVector x;
            double rhs_norm;
            DVector r;
            double r_norm;
            double resid;
        };

        virtual Data init_method(const AnyMatrix& matrix, const DVector& rhs, const DVector& x0);

        void init_preconditioner(const AnyMatrix& matrix) const;
        DVector apply_precondition(const DVector& vector) const;

        bool check_convergence_criterion(const DVector& residual, double rhs_norm, size_t iteration) const;

        size_t get_max_iteration_count(size_t system_size) const;
        void collect_solution_stats(const SolutionStats& stats);

    protected:
        const IterativeSolvingParams m_iterative_solving_params;

    private:
        SolutionStats m_solution_stats;
        const std::unique_ptr<IPreconditioner> m_preconditioner;
    };

} // namespace math::linal