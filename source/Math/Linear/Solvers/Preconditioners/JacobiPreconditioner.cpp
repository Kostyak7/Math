#include <Utils/Math/Linear/Solvers/Preconditioners/JacobiPreconditioner.h>

#include <Utils/Math/Common.h>

#include <stdexcept>

class math::linal::JacobiPreconditioner::Impl : public math::linal::IPreconditioner::IImpl {
public:
    DVector apply(const DVector& x) const override {
        DVector result(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            result[i] = x[i] * m_inverse_diagonal[i];
        }
        return result;
    }

protected:
    DVector m_inverse_diagonal;
};

namespace {

    class BandImpl : public math::linal::JacobiPreconditioner::Impl {
    public:
        BandImpl(const math::linal::BandMatrix& matrix) {
            const size_t n = matrix.get_height();
            m_inverse_diagonal.resize(n);

            for (size_t i = 0; i < n; ++i) {
                const double diag = matrix[i][0];

                if (math::dcmp(diag) == 0) {
                    throw std::runtime_error("Zero diagonal element found in matrix for Jacobi preconditioner");
                }

                m_inverse_diagonal[i] = 1.0 / diag;
            }
        }
    };

    class DenseImpl : public math::linal::JacobiPreconditioner::Impl {
    public:
        DenseImpl(const math::linal::DenseMatrix& matrix) {
            const size_t n = matrix.get_height();
            m_inverse_diagonal.resize(n);

            for (size_t i = 0; i < n; ++i) {
                const double diag = matrix[i][i];

                if (::math::dcmp(diag) == 0) {
                    throw std::runtime_error("Zero diagonal element found in matrix for Jacobi preconditioner");
                }

                m_inverse_diagonal[i] = 1.0 / diag;
            }
        }
    };

} // namespace

void math::linal::JacobiPreconditioner::init(const AnyMatrixConstRef& matrix) {
    if (std::holds_alternative<std::reference_wrapper<const BandMatrix>>(matrix)) {
        m_impl = std::make_unique<BandImpl>(std::get<std::reference_wrapper<const BandMatrix>>(matrix).get());
    }
    else if (std::holds_alternative<std::reference_wrapper<const DenseMatrix>>(matrix)) {
        m_impl = std::make_unique<DenseImpl>(std::get<std::reference_wrapper<const DenseMatrix>>(matrix).get());
    }
    else {
        throw std::runtime_error("Unsupported matrix type for Jacobi preconditioner");
    }
}
