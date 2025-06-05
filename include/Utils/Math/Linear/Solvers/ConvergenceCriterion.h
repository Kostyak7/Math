#pragma once 

#include "Interfaces/IConvergenceCriterion.h"

#include <Utils/Math/Common.h>

namespace math::linal {

    class MATH_EXPORT RelativeResidualCriterion : public IConvergenceCriterion {
    public:
        explicit RelativeResidualCriterion(double tolerance = TOLERANCE)
            : m_tolerance(tolerance) 
        {
        }

        bool is_converged(const DVector& residual, double rhs_norm, size_t /*iteration*/) const override {
            if (rhs_norm < m_tolerance)
                rhs_norm = 1.0;
            return residual.norm() / rhs_norm < m_tolerance;
        }

    private:
        double m_tolerance;
    };

    class AbsoluteResidualCriterion : public IConvergenceCriterion {
    public:
        explicit AbsoluteResidualCriterion(double tolerance) 
            : m_tolerance(tolerance)
        {
        }

        bool is_converged(const DVector& residual, double /*rhs_norm*/, size_t /*iteration*/) const override {
            return residual.norm() < m_tolerance;
        }

    private:
        double m_tolerance;
    };

    class MaxComponentCriterion : public IConvergenceCriterion {
    public:
        explicit MaxComponentCriterion(double tolerance) 
            : m_tolerance(tolerance) 
        {
        }

        bool is_converged(const DVector& residual, double /*rhs_norm*/, size_t /*iteration*/) const override {
            double max_component = 0.0;
            for (double val : residual)
                max_component = std::max(max_component, std::abs(val));
            return max_component < m_tolerance;
        }

    private:
        double m_tolerance;
    };

    class CompositeCriterion : public IConvergenceCriterion {
    public:
        CompositeCriterion(std::vector<std::shared_ptr<IConvergenceCriterion>> criterions)
            : m_criterions(std::move(criterions))
        {
        }

        bool is_converged(const DVector& residual, double rhs_norm, size_t iteration) const override {
            for (const auto& criterion : m_criterions) {
                if (criterion->is_converged(residual, rhs_norm, iteration))
                    return true;
            }
            return false;
        }

    private:
        std::vector<std::shared_ptr<IConvergenceCriterion>> m_criterions;
    };

} // namespace math::linal