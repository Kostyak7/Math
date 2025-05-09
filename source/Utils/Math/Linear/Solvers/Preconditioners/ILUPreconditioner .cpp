#include "ILUPreconditioner.h"

#include <Utils/Math/Common.h>

#include <stdexcept>

class math::linal::ILUPreconditioner::Impl : public math::linal::IPreconditioner::IImpl {
public:
    FVector apply(const FVector& x) const override {
        // ...
        return x;
    }

protected:

};

namespace {

    class BandImpl : public math::linal::ILUPreconditioner::Impl {
        using Matrix = math::linal::BandMatrix;
        using Vector = math::linal::FVector;

    public:
        BandImpl(const Matrix& matrix) {
            const size_t n = matrix.get_width();
            isl = matrix.get_isl();

            L.resize(n);
            U.resize(n);

            for (size_t i = 0; i < n; ++i) {
                // U 
                size_t u_size = std::min(isl, n - i);
                U[i].resize(u_size);
                for (size_t j = 0; j < u_size; ++j) {
                    U[i][j] = matrix[i][j];
                }

                // L 
                size_t l_size = i < isl ? i : isl;
                L[i].resize(l_size, 0.0);
                if (i > 0 && L[i].size() > 0) {
                    L[i][L[i].size() - 1] = 1.0; 
                }
            }

            // ILU-разложение 
            for (size_t k = 0; k < n; ++k) {
                double& u_kk = U[k][0];
                double sum = 0.0;

                size_t start_m = k > isl ? k - isl : 0;
                for (size_t m = start_m; m < k; ++m) {
                    size_t l_pos = m - (k - L[k].size());
                    size_t u_pos = k - m;

                    if (l_pos < L[k].size() && u_pos < U[m].size()) {
                        sum += L[k][l_pos] * U[m][u_pos];
                    }
                }
                u_kk -= sum;

                if (math::dcmp(u_kk) == 0) 
                    u_kk = math::TOLERANCE;

                // Обновляем элементы U в строке k
                size_t max_j = std::min(k + isl, n);
                for (size_t j = k + 1; j < max_j; ++j) {
                    double& u_kj = U[k][j - k];
                    sum = 0.0;

                    start_m = k > isl ? k - isl : 0;
                    for (size_t m = start_m; m < k; ++m) {
                        size_t l_pos = m - (k - L[k].size());
                        size_t u_pos = j - m;

                        if (l_pos < L[k].size() && u_pos < U[m].size()) {
                            sum += L[k][l_pos] * U[m][u_pos];
                        }
                    }
                    u_kj -= sum;
                }

                // Обновляем элементы L в столбце k
                size_t max_i = std::min(k + isl, n);
                for (size_t i = k + 1; i < max_i; ++i) {
                    size_t l_pos = k - (i - L[i].size());
                    if (l_pos >= L[i].size()) continue;

                    double& l_ik = L[i][l_pos];
                    sum = 0.0;

                    start_m = i > isl ? i - isl : 0;
                    start_m = std::max(start_m, k > isl ? k - isl : 0);
                    for (size_t m = start_m; m < k; ++m) {
                        size_t l_i_pos = m - (i - L[i].size());
                        size_t u_pos = k - m;

                        if (l_i_pos < L[i].size() && u_pos < U[m].size()) {
                            sum += L[i][l_i_pos] * U[m][u_pos];
                        }
                    }

                    if (U[k].size() > 0) {
                        l_ik = (matrix.get(i, k) - sum) / u_kk;
                    }
                }
            }
        }

        Vector apply(const Vector& x) const override {
            const size_t n = x.size();
            Vector y(n);
            Vector z(n);

            // Ly = x 
            for (size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                size_t start_j = i > isl ? i - isl : 0;
                for (size_t j = start_j; j < i; ++j) {
                    size_t pos = j - (i - L[i].size());
                    if (pos < L[i].size()) {
                        sum += L[i][pos] * y[j];
                    }
                }
                y[i] = x[i] - sum;
            }

            // Uz = y 
            for (int i = n - 1; i >= 0; --i) {
                double sum = 0.0;
                size_t end_j = std::min(i + isl, n);
                for (size_t j = i + 1; j < end_j; ++j) {
                    size_t pos = j - i;
                    if (pos < U[i].size()) {
                        sum += U[i][pos] * z[j];
                    }
                }
                if (U[i].size() > 0) {
                    z[i] = (y[i] - sum) / U[i][0];
                }
            }

            return z;
        }

    private:
        std::vector<Vector> L; // Нижняя треугольная матрица (включая диагональ)
        std::vector<Vector> U; // Верхняя треугольная матрица (включая диагональ)
        size_t isl;
    };

    class DenseImpl : public math::linal::ILUPreconditioner::Impl {
    public:
        DenseImpl(const math::linal::DenseMatrix& matrix) {
            // ...
        }

    private:

    };

    class SparseImpl : public math::linal::ILUPreconditioner::Impl {
    public:
        SparseImpl(const math::linal::SparseMatrix& matrix) {
            // ...
        }

    private:

    };

} // namespace

void math::linal::ILUPreconditioner::init(const AnyMatrix& matrix) {
    if (std::holds_alternative<BandMatrix>(matrix)) {
        m_impl = std::make_unique<BandImpl>(std::get<BandMatrix>(matrix));
    } 
    else if (std::holds_alternative<DenseMatrix>(matrix)) {
        m_impl = std::make_unique<DenseImpl>(std::get<DenseMatrix>(matrix));
    }
    else if (std::holds_alternative<SparseMatrix>(matrix)) {
        m_impl = std::make_unique<SparseImpl>(std::get<SparseMatrix>(matrix));
    }    
    else {
        throw std::runtime_error("Unsupported matrix type for ILU preconditioner");
    }
}