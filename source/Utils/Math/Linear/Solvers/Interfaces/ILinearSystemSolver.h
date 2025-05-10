#pragma once 

#include <Utils/Math/Linear/AnyMatrix.h>

namespace math::linal {

    class ILinearSystemSolver {
    public:
        struct Params {
            bool throw_exceptions = true;
        };

        ILinearSystemSolver(const Params& params = {});

        virtual DVector solve(const AnyMatrix& matrix, const DVector& rhs, const DVector& x0 = {}) = 0;
        virtual ~ILinearSystemSolver() = default;

        Params get_params() const;

    protected:
        virtual bool slae_check(const AnyMatrix& matrix, const DVector& rhs) const;

    protected:
        Params m_params;
    };

} // namespace math::linal