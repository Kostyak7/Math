#include "DenseMatrix.h"

#include <Utils/Math/Common.h>

#include <stdexcept>

math::linal::DenseMatrix::DenseMatrix(size_t height, size_t width, const value_type& default_value)
    : vector(height, DVector(width, default_value))
{
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator*=(value_type scalar) {
    for (auto& vec : *this) {
        vec *= scalar;
    }
    return *this;
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator/=(value_type scalar) {
    return *this *= (1 / scalar);
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator+=(const DenseMatrix& other) {
    if (get_width() != other.get_width() || get_height() != other.get_height())
        throw std::invalid_argument("The matrices must be the same size");
    for (size_t r = 0; r < get_height(); ++r) {
        (*this)[r] += other[r];
    }
    return *this;
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator-=(const DenseMatrix& other) {
    if (get_width() != other.get_width() || get_height() != other.get_height())
        throw std::invalid_argument("The matrices must be the same size");
    for (size_t r = 0; r < get_height(); ++r) {
        (*this)[r] -= other[r];
    }
    return *this;
}

bool math::linal::DenseMatrix::is_positive_definite() const {
    if (!is_square()) 
        return false;
    if (!is_symmetrical())
        return IMatrix::is_positive_definite();

    const size_t n = get_width();
    DenseMatrix L(n, n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            value_type sum = 0;
            for (size_t k = 0; k < j; ++k) {
                sum += L[i][k] * L[j][k];
            }

            if (i == j) {
                value_type diag = (*this)[i][i] - sum;
                if (diag <= 0) 
                    return false; 

                L[i][i] = std::sqrt(diag);
            }
            else {
                L[i][j] = ((*this)[i][j] - sum) / L[j][j];
            }
        }
    }
    return true;
}

math::linal::DenseMatrix::value_type math::linal::DenseMatrix::get(size_t row, size_t col) const {
    if (row >= get_height() || col >= get_width())
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    return (*this)[row][col];
}

void math::linal::DenseMatrix::set(size_t row, size_t col, value_type value) {
    if (row >= get_height() || col >= get_width())
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    (*this)[row][col] = value;
}

size_t math::linal::DenseMatrix::get_width() const {
    if (empty()) {
        return 0;
    }
    return at(0).size();
}

size_t math::linal::DenseMatrix::get_height() const {
    return size();
}

math::linal::DenseMatrix::value_type math::linal::DenseMatrix::det() const {
    if (!is_square()) {
        throw std::runtime_error("Determinant is defined only for square matrices");
    }

    size_t n = get_width();
    if (n == 0) return 0.0;
    if (n == 1) return (*this)[0][0];
    if (n == 2) return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];

    DenseMatrix temp(*this);
    value_type det = 1.0;

    for (size_t i = 0; i < n; ++i) {
        size_t pivot = i;
        for (size_t j = i + 1; j < n; ++j) {
            if (fabs(temp[j][i]) > fabs(temp[pivot][i])) {
                pivot = j;
            }
        }

        if (pivot != i) {
            std::swap(temp[i], temp[pivot]);
            det = -det;
        }

        if (temp[i][i] == 0) 
            return 0.0;

        det *= temp[i][i];

        for (size_t j = i + 1; j < n; ++j) {
            value_type factor = temp[j][i] / temp[i][i];
            for (size_t k = i + 1; k < n; ++k) {
                temp[j][k] -= factor * temp[i][k];
            }
        }
    }

    return det;
}

std::vector<std::pair<size_t, math::linal::DenseMatrix::complex_value_type>> math::linal::DenseMatrix::get_eigenvalues() const {
    std::vector<std::pair<size_t, complex_value_type>> eigenvalues;
    // ...
    return eigenvalues;
}

std::vector<std::pair<math::linal::DenseMatrix::complex_value_type, std::vector<math::linal::DVector>>> math::linal::DenseMatrix::get_eigenvectors() const {
    std::vector<std::pair<complex_value_type, std::vector<math::linal::DVector>>> eigenvectors;
    auto eigenvalues = get_eigenvalues();
    // ...
    return eigenvectors;
}

math::linal::DenseMatrix math::linal::DenseMatrix::identity_matrix(size_t n) {
    DenseMatrix res(n, n);
    for (size_t i = 0; i < n; ++i) {
        res[i][i] = 1.0;
    }
    return res;
}

math::linal::DenseMatrix math::linal::DenseMatrix::elementary_matrix_unit(size_t n, size_t m, size_t i, size_t j, double value) {
    DenseMatrix res(n, m);
    res.set(i, j, value);
    return res;
}

math::linal::DenseMatrix math::linal::operator*(const DenseMatrix& matrix, DenseMatrix::value_type scalar) {
    DenseMatrix res(matrix);
    res *= scalar;
    return res;
}

math::linal::DenseMatrix math::linal::operator*(DenseMatrix::value_type scalar, const DenseMatrix& matrix) {
    return matrix * scalar;
}

math::linal::DenseMatrix math::linal::operator/(const DenseMatrix& matrix, DenseMatrix::value_type scalar) {
    return matrix * (1. / scalar);
}

bool math::linal::operator==(const DenseMatrix& m1, const DenseMatrix& m2) {
    if (&m1 == &m2)
        return true;
    if (m1.get_width() != m2.get_width() || m1.get_height() != m2.get_height())
        return false;
    for (size_t r = 0; r < m1.get_height(); ++r) {
        if (m1[r] != m2[r]) {
            return false;
        }
    }
    return true;
}

bool math::linal::operator!=(const DenseMatrix& m1, const DenseMatrix& m2) {
    return !(m1 == m2);
}

math::linal::DVector math::linal::operator*(const DenseMatrix& matrix, const DVector& vector) {
    if (matrix.get_width() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }
    const size_t n = matrix.get_height();
    DVector res(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        res[i] = matrix[i].dot(vector);
    }
    return res;
}

math::linal::DVector math::linal::operator*(const DVector& vector, const DenseMatrix& matrix) {
    if (matrix.get_height() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }
    const size_t n = matrix.get_width();
    DVector res(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < matrix.get_height(); ++j) {
            res[i] += vector[j] * matrix[j][i];
        }
    }
    return res;
}

math::linal::DenseMatrix math::linal::operator*(const DenseMatrix& m1, const DenseMatrix& m2) {
    if (m1.get_width() != m2.get_height()) {
        throw std::invalid_argument("The matrices must be the same size");
    }
    DenseMatrix res(m1.get_height(), DVector(m2.get_width()));
    for (size_t r = 0; r < res.get_height(); ++r) {
        for (size_t c = 0; c < res.get_width(); ++c) {
            DenseMatrix::value_type s = 0.0;
            for (size_t i = 0; i < m1.get_width(); ++i) {
                s += m1[r][i] * m2[i][c];
            }
            res[r][c] = s;
        }
    }
    return res;
}

math::linal::DenseMatrix math::linal::operator+(const DenseMatrix& m1, const DenseMatrix& m2) {
    if (m1.get_width() != m2.get_width() || m1.get_height() != m2.get_height())
        throw std::invalid_argument("The matrices must be the same size");
    DenseMatrix res(m1);
    res += m2;
    return res;
}

math::linal::DenseMatrix math::linal::operator-(const DenseMatrix& m1, const DenseMatrix& m2) {
    if (m1.get_width() != m2.get_width() || m1.get_height() != m2.get_height())
        throw std::invalid_argument("The matrices must be the same size");
    DenseMatrix res(m1);
    res -= m2;
    return res;
}

math::linal::DenseMatrix math::linal::operator-(const DenseMatrix& matrix) {
    return -1 * matrix;
}

math::linal::DenseMatrix math::linal::inversed(const DenseMatrix& matrix) {
    if (matrix.is_empty()) 
        throw std::invalid_argument("Matrix is empty");
    if (!matrix.is_sqaure()) {
        throw std::runtime_error("Inverse is defined only for square matrices");
    }

    const size_t n = matrix.get_width();
    DenseMatrix inverse = DenseMatrix::identity_matrix(n);
    DenseMatrix temp(matrix);

    // Прямой ход метода Гаусса
    for (size_t i = 0; i < n; ++i) {
        size_t pivot = i;
        for (size_t j = i + 1; j < n; ++j) {
            if (std::abs(temp[j][i]) > std::abs(temp[pivot][i])) {
                pivot = j;
            }
        }

        if (pivot != i) {
            std::swap(temp[i], temp[pivot]);
            std::swap(inverse[i], inverse[pivot]);
        }

        if (temp[i][i] == 0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted");
        }

        DenseMatrix::value_type diag = temp[i][i];
        for (size_t j = 0; j < n; ++j) {
            temp[i][j] /= diag;
            inverse[i][j] /= diag;
        }

        for (size_t k = 0; k < n; ++k) {
            if (k != i && temp[k][i] != 0) {
                DenseMatrix::value_type factor = temp[k][i];
                for (size_t j = 0; j < n; ++j) {
                    temp[k][j] -= factor * temp[i][j];
                    inverse[k][j] -= factor * inverse[i][j];
                }
            }
        }
    }

    return inverse;
}

math::linal::DenseMatrix math::linal::transposed(const DenseMatrix& matrix) {
    DenseMatrix res(matrix.get_width(), matrix.get_height());
    for (size_t r = 0; r < matrix.get_height(); ++r) {
        for (size_t c = 0; c < matrix.get_width(); ++c) {
            res[c][r] = matrix[r][c];
        }
    }
    return res;
}
