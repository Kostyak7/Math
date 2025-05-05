#include "GaussSystemSolver.h"

#include <Utils/Math/Common.h>

#include <stdexcept>

namespace {

    void forward_elimination(math::linal::BandMatrix& matrix, math::linal::FVector& rhs) {
        int s, m, k, j;
        const size_t n = rhs.size();
        const size_t isl = matrix.get_isl();
        double zn;
        thread_local std::vector<double> current_row;
        current_row.resize(isl, 0.0);

        for (int r = 0; r < n; ++r) {
            rhs[r] /= matrix[r][0];
            if (r == n - 1)
                break;
            zn = matrix[r][0];
            for (s = 1; s < isl; ++s) {
                current_row[s] = matrix[r][s];
                if (math::dcmp(current_row[s]) == 0)
                    continue;
                matrix[r][s] = current_row[s] / zn;
            }
            for (m = 1; m < isl; ++m) {
                zn = current_row[m];
                if (math::dcmp(zn) == 0)
                    continue;
                k = r + m;
                if (k > n - 1)
                    continue;
                j = 0;
                for (s = m; s < isl; ++s, ++j) {
                    matrix[k][j] -= zn * matrix[r][s];
                    if (math::dcmp(matrix[k][0]) == 0)
                        matrix[k][0] = math::TOLERANCE; 
                }
                rhs[k] -= zn * rhs[r];
            }
        }
        return;
    }

    void backward_substitution(const math::linal::BandMatrix& matrix, math::linal::FVector& rhs) {
        const size_t n = rhs.size();
        const size_t isl = matrix.get_isl();

        for (int r = n - 2; r >= 0; --r) {
            for (int s = 1; s < isl; ++s) {
                int m = r + s;
                if (m > n - 1)
                    continue;
                rhs[r] -= matrix[r][s] * rhs[m];
            }
        }
    }

    void forward_elimination(math::linal::DenseMatrix& matrix, math::linal::FVector& rhs) {
        // ...
    }

    void backward_substitution(const math::linal::DenseMatrix& matrix, math::linal::FVector& rhs) {
        // ...
    }

    void forward_elimination(math::linal::SparseMatrix& matrix, math::linal::FVector& rhs) {
        // ...
    }

    void backward_substitution(const math::linal::SparseMatrix& matrix, math::linal::FVector& rhs) {
        // ...
    }

} // namespace 

math::linal::GaussLinearSystemSolver::GaussLinearSystemSolver(const Params& params) 
    : ILinearSystemSolver(params)
{
}

math::linal::FVector math::linal::GaussLinearSystemSolver::solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& /*x0*/) {
    return std::visit([this, &rhs](const auto& matrix) -> FVector {
        if (!slae_check(matrix, rhs)) {
            return {};
        }

        auto augmented_matrix = matrix;
        FVector solution = rhs;

        forward_elimination(augmented_matrix, solution);

        backward_substitution(augmented_matrix, solution);

        return solution;
        
        }, matrix);
}
