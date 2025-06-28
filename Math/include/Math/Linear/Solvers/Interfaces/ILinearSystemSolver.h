#pragma once 

#include <Math/Linear/AnyMatrix.h>

namespace math::linal {

    class MATH_EXPORT ILinearSystemSolver {
    public:
        struct Params {
            bool throw_exceptions = true;
        };

        ILinearSystemSolver(const Params& params = {});

        virtual DVector solve(const AnyMatrixConstRef& matrix, const DVector& rhs, const DVector& x0 = {}) = 0;
        virtual ~ILinearSystemSolver() = default;

        Params get_params() const;

    protected:
        virtual bool slae_check(const IMatrix& matrix, const DVector& rhs) const;

    protected:
        Params m_params;
    };

} // namespace math::linal