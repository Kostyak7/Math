#pragma once 

#include "BandMatrix.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"

#include <variant>

namespace math::linal {

    using AnyMatrix = std::variant<DenseMatrix, BandMatrix, SparseMatrix>;

    AnyMatrix operator*(const AnyMatrix& matrix, IMatrix::value_type scalar);
    AnyMatrix operator*(IMatrix::value_type scalar, const AnyMatrix& matrix);
    AnyMatrix operator/(const AnyMatrix& matrix, IMatrix::value_type scalar);

    DVector operator*(const AnyMatrix& matrix, const DVector& vector);
    DVector operator*(const DVector& vector, const AnyMatrix& matrix);

    AnyMatrix operator*(const AnyMatrix& m1, const AnyMatrix& m2);

    AnyMatrix operator+(const AnyMatrix& m1, const AnyMatrix& m2);
    AnyMatrix operator-(const AnyMatrix& m1, const AnyMatrix& m2);

} // namespace math::linal

