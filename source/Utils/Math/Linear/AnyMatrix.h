#pragma once 

#include "BandMatrix.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"

#include <variant>
#include <functional>

namespace math::linal {

    using AnyMatrix = std::variant<
        BandMatrix, 
        DenseMatrix, 
        SparseMatrix
    >;
    using AnyMatrixRef = std::variant<
        std::reference_wrapper<BandMatrix>,
        std::reference_wrapper<DenseMatrix>,        
        std::reference_wrapper<SparseMatrix>
    >;
    using AnyMatrixConstRef = std::variant<
        std::reference_wrapper<const BandMatrix>,
        std::reference_wrapper<const DenseMatrix>,
        std::reference_wrapper<const SparseMatrix>
    >;

    AnyMatrix operator*(const AnyMatrix& matrix, IMatrix::value_type scalar);
    AnyMatrix operator*(IMatrix::value_type scalar, const AnyMatrix& matrix);
    AnyMatrix operator/(const AnyMatrix& matrix, IMatrix::value_type scalar);

    DVector operator*(const AnyMatrix& matrix, const DVector& vector);
    DVector operator*(const DVector& vector, const AnyMatrix& matrix);

    AnyMatrix operator*(const AnyMatrix& m1, const AnyMatrix& m2);

    AnyMatrix operator+(const AnyMatrix& m1, const AnyMatrix& m2);
    AnyMatrix operator-(const AnyMatrix& m1, const AnyMatrix& m2);

} // namespace math::linal

