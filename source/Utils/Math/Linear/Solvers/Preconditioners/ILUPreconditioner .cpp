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
        BandImpl(const Matrix& matrix) 
            : m_LU(matrix)
        {
            const size_t n = m_LU.get_width();
            const size_t isl = m_LU.get_isl();

            for (size_t i = 0; i < n; ++i) {
                size_t j_min = (i >= isl - 1) ? i - (isl - 1) : 0;

                for (size_t j = j_min; j < i; ++j) {
                    size_t j_offset = i - j;
                    if (j_offset >= isl) continue; // вне ленты

                    double sum = m_LU[i][j_offset];

                    // k от j_min до j-1
                    for (size_t k = j_min; k < j; ++k) {
                        size_t ik = i - k;
                        size_t jk = j - k;
                        if (ik >= isl || jk >= isl) continue;

                        sum -= m_LU[i][ik] * m_LU[k][jk];
                    }

                    // l_ij = sum / u_jj
                    m_LU[i][j_offset] = sum / m_LU[j][0];
                }

                // диагональ u_ii
                double sum = m_LU[i][0];
                for (size_t k = j_min; k < i; ++k) {
                    size_t ik = i - k;
                    if (ik >= isl) continue;

                    sum -= m_LU[i][ik] * m_LU[i][ik]; // l_ik * u_ik
                }

                m_LU[i][0] = sum;
            }
        }

        Vector apply(const Vector& x) const override {
            const size_t n = m_LU.get_width();
            const size_t isl = m_LU.get_isl();

            Vector y(n, 0.0);
            Vector z(n, 0.0);

            // Прямой ход: L * y = x
            for (size_t i = 0; i < n; ++i) {
                double sum = x[i];
                size_t j_min = (i >= isl - 1) ? i - (isl - 1) : 0;

                for (size_t j = j_min; j < i; ++j) {
                    size_t offset = i - j;
                    if (offset >= isl) continue;
                    sum -= m_LU[i][offset] * y[j];
                }

                y[i] = sum;
            }

            // Обратный ход: U * z = y
            for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
                double sum = y[i];

                size_t j_max = std::min(n, i + isl);
                for (size_t j = i + 1; j < j_max; ++j) {
                    size_t offset = j - i;
                    if (offset >= isl) continue;
                    sum -= m_LU[i][offset] * z[j];
                }

                z[i] = sum / m_LU[i][0];
            }

            return z;
        }

    private:
        Matrix m_LU;
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