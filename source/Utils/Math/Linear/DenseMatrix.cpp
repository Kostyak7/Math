#include "DenseMatrix.h"

#include <Utils/Math/Common.h>

#include <stdexcept>

math::linal::DenseMatrix::ProxyVector::ProxyVector(FVector& vector)
    : m_vector(vector)
{
}

const math::linal::DenseMatrix::value_type& math::linal::DenseMatrix::ProxyVector::operator[](size_t i) const {
    return m_vector[i];
}

math::linal::DenseMatrix::value_type& math::linal::DenseMatrix::ProxyVector::operator[](size_t i) {
    return m_vector[i];
}

math::linal::DenseMatrix::DenseMatrix(size_t height, size_t width, const value_type& default_value)
    : m_data(height, FVector(width, default_value))
{
}

math::linal::DenseMatrix::DenseMatrix(size_t height, const FVector& default_value)
    : m_data(height, default_value)
{
}

math::linal::DenseMatrix::DenseMatrix(const DenseMatrix& matrix)
    : m_data(matrix.m_data)
{
}

math::linal::DenseMatrix::DenseMatrix(DenseMatrix&& matrix) noexcept {
    swap(std::move(matrix));
}

math::linal::DenseMatrix::DenseMatrix(std::initializer_list<value_type> list)
    : m_data(list.begin(), list.end())
{
}

math::linal::DenseMatrix::DenseMatrix(std::initializer_list<std::initializer_list<value_type>> list) {
    for (const auto& row : list) {
        m_data.insert(m_data.end(), row.begin(), row.end());
    }
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator=(const DenseMatrix& other) {
    if (this == &other)
        return *this;
    m_data = other.m_data;
    return *this;
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator=(DenseMatrix&& other) {
    if (this == &other)
        return *this;
    swap(other);
    return *this;

}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator*=(value_type scalar) noexcept {
    for (auto& vec : m_data) {
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
        m_data[r] += other[r];
    }
    return *this;
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator-=(const DenseMatrix& other) {
    if (get_width() != other.get_width() || get_height() != other.get_height())
        throw std::invalid_argument("The matrices must be the same size");
    for (size_t r = 0; r < get_height(); ++r) {
        m_data[r] -= other[r];
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
                sum += L.m_data[i][k] * L.m_data[j][k];
            }

            if (i == j) {
                value_type diag = m_data[i][i] - sum;
                if (diag <= 0) 
                    return false; 

                L.m_data[i][i] = std::sqrt(diag);
            }
            else {
                L.m_data[i][j] = (m_data[i][j] - sum) / L.m_data[j][j];
            }
        }
    }
    return true;
}

math::linal::DenseMatrix::value_type math::linal::DenseMatrix::get(size_t row, size_t col) const {
    if (row >= get_height() || col >= get_width())
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    return m_data[row][col];
}

void math::linal::DenseMatrix::set(size_t row, size_t col, value_type value) {
    if (row >= get_height() || col >= get_width())
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    m_data[row][col] = value;
}

size_t math::linal::DenseMatrix::get_width() const {
    return m_data.at(0).size();
}

size_t math::linal::DenseMatrix::get_height() const {
    return m_data.size();
}

const math::linal::FVector& math::linal::DenseMatrix::operator[](size_t row) const {
    return m_data[row];
}

math::linal::DenseMatrix::ProxyVector math::linal::DenseMatrix::operator[](size_t row) {
    return { m_data[row] };
}

void math::linal::DenseMatrix::clear() {
    m_data.clear();
}

math::linal::DenseMatrix::value_type math::linal::DenseMatrix::det() const {
    if (!is_square()) {
        throw std::runtime_error("Determinant is defined only for square matrices");
    }

    size_t n = get_width();
    if (n == 0) return 0.0;
    if (n == 1) return m_data[0][0];
    if (n == 2) return m_data[0][0] * m_data[1][1] - m_data[0][1] * m_data[1][0];

    DenseMatrix temp(*this);
    value_type det = 1.0;

    for (size_t i = 0; i < n; ++i) {
        size_t pivot = i;
        for (size_t j = i + 1; j < n; ++j) {
            if (fabs(temp.m_data[j][i]) > fabs(temp.m_data[pivot][i])) {
                pivot = j;
            }
        }

        if (pivot != i) {
            std::swap(temp.m_data[i], temp.m_data[pivot]);
            det = -det;
        }

        if (temp.m_data[i][i] == 0) 
            return 0.0;

        det *= temp.m_data[i][i];

        for (size_t j = i + 1; j < n; ++j) {
            value_type factor = temp.m_data[j][i] / temp.m_data[i][i];
            for (size_t k = i + 1; k < n; ++k) {
                temp.m_data[j][k] -= factor * temp.m_data[i][k];
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

std::vector<std::pair<math::linal::DenseMatrix::complex_value_type, std::vector<math::linal::FVector>>> math::linal::DenseMatrix::get_eigenvectors() const {
    std::vector<std::pair<complex_value_type, std::vector<math::linal::FVector>>> eigenvectors;
    auto eigenvalues = get_eigenvalues();
    // ...
    return eigenvectors;
}

math::linal::DenseMatrix math::linal::DenseMatrix::identity_matrix(size_t n) {
    DenseMatrix res(n, n);
    for (size_t i = 0; i < n; ++i) {
        res.m_data[i][i] = 1.0;
    }
    return res;
}

const std::vector<math::linal::FVector>& math::linal::DenseMatrix::data() const {
    return m_data;
}

std::vector<math::linal::FVector>& math::linal::DenseMatrix::data() {
    return m_data;
}

void math::linal::DenseMatrix::swap(DenseMatrix& other) noexcept {
    std::swap(m_data, other.m_data);
}

void math::linal::swap(DenseMatrix& m1, DenseMatrix& m2) noexcept {
    m1.swap(m2);
}

math::linal::DenseMatrix math::linal::operator*(const DenseMatrix& matrix, DenseMatrix::value_type scalar) {
    DenseMatrix res(matrix);
    res *= scalar;
    return res;
}

math::linal::DenseMatrix math::linal::operator*(DenseMatrix::value_type scalar, const DenseMatrix& matrix) {
    DenseMatrix res(matrix);
    res *= scalar;
    return res;
}

math::linal::DenseMatrix math::linal::operator/(const DenseMatrix& matrix, DenseMatrix::value_type scalar) {
    DenseMatrix res(matrix);
    res /= scalar;
    return res;
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

math::linal::FVector math::linal::operator*(const DenseMatrix& matrix, const FVector& vector) {
    if (matrix.get_width() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }
    const size_t n = matrix.get_height();
    FVector res(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        res[i] = matrix[i].dot(vector);
    }
    return res;
}

math::linal::FVector math::linal::operator*(const FVector& vector, const DenseMatrix& matrix) {
    if (matrix.get_height() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }
    const size_t n = matrix.get_width();
    FVector res(n, 0.0);
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
    DenseMatrix res(m1.get_height(), FVector(m2.get_width()));
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
    if (matrix.get_width() != matrix.get_height()) {
        throw std::runtime_error("Inverse is defined only for square matrices");
    }

    const size_t n = matrix.get_width();
    DenseMatrix inverse = DenseMatrix::identity_matrix(n);
    DenseMatrix temp(matrix);

    // Прямой ход метода Гаусса
    for (size_t i = 0; i < n; ++i) {
        size_t pivot = i;
        for (size_t j = i + 1; j < n; ++j) {
            if (std::abs(temp.m_data[j][i]) > std::abs(temp.m_data[pivot][i])) {
                pivot = j;
            }
        }

        if (pivot != i) {
            std::swap(temp.m_data[i], temp.m_data[pivot]);
            std::swap(inverse.m_data[i], inverse.m_data[pivot]);
        }

        if (temp.m_data[i][i] == 0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted");
        }

        DenseMatrix::value_type diag = temp.m_data[i][i];
        for (size_t j = 0; j < n; ++j) {
            temp.m_data[i][j] /= diag;
            inverse.m_data[i][j] /= diag;
        }

        // Исключение элементов в столбце
        for (size_t k = 0; k < n; ++k) {
            if (k != i && temp.m_data[k][i] != 0) {
                DenseMatrix::value_type factor = temp.m_data[k][i];
                for (size_t j = 0; j < n; ++j) {
                    temp.m_data[k][j] -= factor * temp.m_data[i][j];
                    inverse.m_data[k][j] -= factor * inverse.m_data[i][j];
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
