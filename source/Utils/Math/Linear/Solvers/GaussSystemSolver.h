#pragma once 

#include "Interfaces/ILinearSystemSolver.h"

namespace math::linal {

    /**
     @brief Решение системы методом гаусса с обходом нулей.
     @details Решает систему уравнений следующего вида matrix*solution=rhs (для матрицы ленточного вида).
    */
    class GaussLinearSystemSolver final : public ILinearSystemSolver {
    public:
        GaussLinearSystemSolver(const Params& params = {});
        FVector solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& x0 = {}) override;
    };

} // namespace math::linal