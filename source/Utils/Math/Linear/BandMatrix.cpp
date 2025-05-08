#include "BandMatrix.h"

#include <Utils/Math/Common.h>

#include <stdexcept>


math::linal::BandMatrix& math::linal::BandMatrix::operator*=(value_type scalar) noexcept {
    for (auto& vec : *this) {
        vec *= scalar;
    }
    return *this;
}

math::linal::BandMatrix& math::linal::BandMatrix::operator/=(value_type scalar) {
    return *this *= (1 / scalar);
}

math::linal::BandMatrix& math::linal::BandMatrix::operator+=(const BandMatrix& other) {
    if (get_height() != other.get_height() || get_width() != other.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }

    if (get_isl() < other.get_isl()) {
        set(0, other.get_isl() - 1, 0.);
    }
    const size_t isl = get_isl();
    for (size_t i = 0; i < get_height(); ++i) {
        for (size_t j = 0; j < isl; ++j) {
            (*this)[i][j] += other[i][j];
        }
    }

    return *this;
}

math::linal::BandMatrix& math::linal::BandMatrix::operator-=(const BandMatrix& other) {
    if (get_height() != other.get_height() || get_width() != other.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }

    if (get_isl() < other.get_isl()) {
        set(0, other.get_isl() - 1, 0.);
    }
    const size_t isl = get_isl();
    for (size_t i = 0; i < get_height(); ++i) {
        for (size_t j = 0; j < isl; ++j) {
            (*this)[i][j] -= other[i][j];
        }
    }
    return *this;
}

bool math::linal::BandMatrix::is_zero() const {
    for (const auto& vec : *this) {
        if (!vec.is_zero()) {
            return false;
        }
    }
    return true;
}

bool math::linal::BandMatrix::is_identity() const {
    if (!is_square())
        return false;
    for (const auto& row : *this) {
        if (dcmp(row.at(0), 1.) != 0) return false;
        for (size_t i = 1; i < row.size(); ++i) {
            if (dcmp(row.at(i)) != 0)
                return false;
        }
    }
    return true;
}

bool math::linal::BandMatrix::is_diagonal() const {
    if (!is_square())
        return false;
    for (const auto& row : *this) {
        if (dcmp(row.at(0)) == 0) 
            return false;
        for (size_t i = 1; row.size(); ++i) {
            if (dcmp(row.at(i)) != 0)
                return false;
        }
    }
    return true;
}

bool math::linal::BandMatrix::is_symmetrical() const {
    return true;
}

bool math::linal::BandMatrix::is_upper_triangular() const {
    return is_square() && is_trapezoidal();
}

bool math::linal::BandMatrix::is_lower_triangular() const {
    return is_diagonal();
}

bool math::linal::BandMatrix::is_trapezoidal() const {
    return is_diagonal();
}

math::linal::BandMatrix::value_type math::linal::BandMatrix::get(size_t row, size_t col) const {
    if (row >= get_height() || col >= get_width()) {
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    }
    if (col >= get_isl() + row || row >= get_isl() + col) {
        return 0.0;
    }
    if (col < row) {
        return (*this)[col][row - col];
    }
    return (*this)[row][col - row];
}

void math::linal::BandMatrix::set(size_t row, size_t col, value_type value) {
    if (row >= get_height() || col >= get_width()) {
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    }
    if (col >= get_isl() + row) {
        set_isl(col - row + 1);
    }
    else if (row >= get_isl() + col) {
        set_isl(row - col + 1);
    }
    if (col < row) {
        (*this)[col][row] = value;
    }
    (*this)[row][col] = value;
}

size_t math::linal::BandMatrix::get_width() const {
    return size();
}

size_t math::linal::BandMatrix::get_height() const {
    return size();
}

size_t math::linal::BandMatrix::get_isl() const {
    if (empty()) {
        return 0;
    }
    return at(0).size();
}

void math::linal::BandMatrix::set_isl(size_t width) {
    if (width == get_isl() || width == 0) return;
    for (auto& vec : (*this)) {
        vec.resize(width);
    }
}

void math::linal::BandMatrix::shrink() {
    const size_t old_isl = get_isl();
    size_t new_isl = 1;
    for (size_t row = 0; row < get_height(); ++row) {
        size_t zero_count = 0;
        for (size_t col = old_isl - 1; col >= 0; --col) {
            if (dcmp((*this)[row][col]) != 0) {
                break;
            }
            ++zero_count;
        }
        if (zero_count == 0)
            return;
        new_isl = std::max(new_isl, old_isl - zero_count);
    }
    for (auto& vec : *this) {
        vec.resize(new_isl);
    }
}

math::linal::BandMatrix::value_type math::linal::BandMatrix::det() const {
    value_type determinant = 1.0;
    if (is_empty())
        return determinant;

    // ...
    return determinant;
}

std::vector<std::pair<size_t, math::linal::BandMatrix::complex_value_type>> math::linal::BandMatrix::get_eigenvalues() const {
    std::vector<std::pair<size_t, complex_value_type>> eigenvalues;
    // ...
    return eigenvalues;
}

std::vector<std::pair<math::linal::BandMatrix::complex_value_type, std::vector<math::linal::FVector>>> math::linal::BandMatrix::get_eigenvectors() const {
    std::vector<std::pair<complex_value_type, std::vector<math::linal::FVector>>> eigenvectors;
    auto eigenvalues = get_eigenvalues();
    // ...
    return eigenvectors;
}

math::linal::BandMatrix math::linal::BandMatrix::identity_matrix(size_t n, size_t isl) {
    BandMatrix res(n, FVector(isl, 0.0));
    for (auto& row : res) {
        row[0] = 1.0;
    }
    return res;
}

math::linal::BandMatrix math::linal::BandMatrix::elementary_matrix_unit(size_t n, size_t i, size_t j) {
    BandMatrix res(n);
    res.set(i, j, 1.0);
    return res;
}

bool math::linal::operator==(const BandMatrix& m1, const BandMatrix& m2) {
    if (&m1 == &m2)
        return true;
    if (m1.get_height() != m2.get_height() || m1.get_width() != m2.get_width())
        return false;
    for (size_t row = 0; row < m1.get_height(); ++row) {
        const size_t min_isl = std::min(m1.get_isl(), m2.get_isl());
        for (size_t col = 0; col < min_isl; ++col) {
            if (dcmp(m1[row][col], m2[row][col]) != 0) {
                return false;
            }
        }
        if (m1.get_isl() != m2.get_isl()) {
            const FVector& vec = (m1.get_isl() < m2.get_isl()) ? m2[row] : m1[row];
            for (size_t i = min_isl; i < vec.size(); ++i) {
                if (dcmp(vec[i]) != 0) {
                    return false;
                }
            }
        }
    }
    return true;
}

bool math::linal::operator!=(const BandMatrix& m1, const BandMatrix& m2) {
    return !(m1 == m2);
}

math::linal::BandMatrix math::linal::operator*(const BandMatrix& matrix, BandMatrix::value_type scalar) {
    BandMatrix res(matrix);
    res *= scalar;
    return res;
}

math::linal::BandMatrix math::linal::operator*(BandMatrix::value_type scalar, const BandMatrix& matrix) {
    BandMatrix res(matrix);
    res *= scalar;
    return res;
}

math::linal::BandMatrix math::linal::operator/(const BandMatrix& matrix, BandMatrix::value_type scalar) {
    BandMatrix res(matrix);
    res /= scalar;
    return res;
}

math::linal::FVector math::linal::operator*(const BandMatrix& matrix, const FVector& vector) {
    if (matrix.get_width() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }

    size_t n = matrix.get_height();
    size_t w = matrix.get_isl();
    FVector res(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = std::max(0, static_cast<int>(i) - static_cast<int>(w)); j <= std::min(n - 1, i + w); ++j) {
            res[i] += matrix.get(i, j) * vector[j];
        }
    }
    return res;
}

math::linal::FVector math::linal::operator*(const FVector& vector, const BandMatrix& matrix) {
    return matrix * vector;
}

math::linal::BandMatrix math::linal::operator*(const BandMatrix& m1, const BandMatrix& m2) {
    if (m1.get_width() != m2.get_height()) {
        throw std::invalid_argument("The matrices must be the same size");
    }

    size_t n = m1.get_height();
    size_t w = m1.get_isl();  // ширина ленты
    BandMatrix result(n, FVector(w + 1, 0.0)); // новая ленточная матрица

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = std::max(0, (int)i - (int)w); j <= std::min(n - 1, i + w); ++j) {
            BandMatrix::value_type sum = 0.0;
            for (size_t k = std::max(0, (int)i - (int)w); k <= std::min(n - 1, i + w); ++k) {
                sum += m1.get(i, k) * m2.get(k, j);
            }
            result.set(i, j, sum);
        }
    }

    result.shrink();
    return result;
}

math::linal::BandMatrix math::linal::operator+(const BandMatrix& m1, const BandMatrix& m2) {
    if (m1.get_height() != m2.get_height() || m1.get_width() != m2.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }
    BandMatrix res(m1);
    res += m2;
    return res;
}

math::linal::BandMatrix math::linal::operator-(const BandMatrix& m1, const BandMatrix& m2) {
    if (m1.get_height() != m2.get_height() || m1.get_width() != m2.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }
    BandMatrix res(m1);
    res -= m2;
    return res;
}

math::linal::BandMatrix math::linal::operator-(const BandMatrix& matrix) {
    return -1 * matrix;
}

math::linal::BandMatrix math::linal::inversed(const BandMatrix& matrix) {
    if (matrix.is_empty()) {
        throw std::invalid_argument("Matrix is empty");
    }

    // Далее реализовано нахождение разложения Холецкого
    size_t n = matrix.get_height();
    BandMatrix L(n, FVector(n, 0.0));
    // ...
    BandMatrix res(n, FVector(n, 0.0));
    // ...
    return res;
}
