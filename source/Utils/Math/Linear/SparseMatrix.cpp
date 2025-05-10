#include "SparseMatrix.h"

#include <Utils/Math/Common.h>

#include <stdexcept>

math::linal::SparseMatrix::SparseMatrix(size_t width, size_t height) 
    : m_width(width)
    , m_height(height)
{
}

math::linal::SparseMatrix::SparseMatrix(const SparseMatrix& matrix)
    : m_width(matrix.m_width)
    , m_height(matrix.m_height)
    , m_data(matrix.m_data)
{
}

math::linal::SparseMatrix::SparseMatrix(SparseMatrix&& matrix) {
    swap(matrix);
}

math::linal::SparseMatrix& math::linal::SparseMatrix::operator=(const SparseMatrix& other) {
    if (this == &other)
        return *this;
    m_data = other.m_data;
    m_width = other.m_width;
    m_height = other.m_height;
    return *this;
}

math::linal::SparseMatrix& math::linal::SparseMatrix::operator=(SparseMatrix&& other) {
    if (this == &other) 
        return *this;
    swap(other);
    return *this;

}

math::linal::SparseMatrix& math::linal::SparseMatrix::operator*=(value_type scalar) noexcept {
    for (auto& [_, val] : m_data) {
        val *= scalar;
    }
    return *this;
}

math::linal::SparseMatrix& math::linal::SparseMatrix::operator/=(value_type scalar) {
    return *this *= (1. / scalar);
}

math::linal::SparseMatrix& math::linal::SparseMatrix::operator+=(const SparseMatrix& other) {
    if (get_height() != other.get_height() || get_width() != other.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }

    for (const auto& [key, value] : other.m_data) {
        if (dcmp(value) != 0) {
            if (auto it = m_data.find(key); it != m_data.end()) {
                it->second += value;
            }
            else {
                m_data.insert({ key, value });
            }

        }
    }
    return *this;
}

math::linal::SparseMatrix& math::linal::SparseMatrix::operator-=(const SparseMatrix& other) {
    if (get_height() != other.get_height() || get_width() != other.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }

    for (const auto& [key, value] : other.m_data) {
        if (dcmp(value) != 0) {
            if (auto it = m_data.find(key); it != m_data.end()) {
                it->second -= value;
            }
            else {
                m_data.insert({ key, -value });
            }

        }
    }
    return *this;
}

bool math::linal::SparseMatrix::is_zero() const {
    for (const auto& [_, val] : m_data) {
        if (dcmp(val) != 0)
            return false;
    }
    return true;
}

bool math::linal::SparseMatrix::is_identity() const {
    if (!is_square())
        return false;
    for (const auto& [key, value] : m_data) {
        if (key.first != key.second) {
            if (dcmp(value) != 0) {
                return false;
            }
        }
        else if (dcmp(value, 1.) != 0){
            return false;
        }
    }
    return true;
}

bool math::linal::SparseMatrix::is_diagonal() const {
    if (!is_square())
        return false;
    for (const auto& [key, value] : m_data) {
        if (key.first != key.second && dcmp(value) != 0) {
            return false;
        }
    }
    return true;
}

bool math::linal::SparseMatrix::is_symmetrical() const {
    if (!is_square())
        return false;
    for (const auto& [key, value] : m_data) {
        if (key.first != key.second) {
            if (auto it = m_data.find({ key.second, key.first }); it != m_data.end() && dcmp(value, it->second) != 0) {
                return false;
            }
            else if(dcmp(value) != 0) {
                return false;
            }
        }
    }
    return true;
}

bool math::linal::SparseMatrix::is_upper_triangular() const {
    return is_square() && is_trapezoidal();
}

bool math::linal::SparseMatrix::is_lower_triangular() const {
    if (!is_square())
        return false;
    for (const auto& [key, value] : m_data) {
        if (key.first < key.second && dcmp(value) != 0) {
            return false;
        }
    }
    return true;
}

bool math::linal::SparseMatrix::is_trapezoidal() const {
    for (const auto& [key, value] : m_data) {
        if (key.first > key.second && dcmp(value) != 0) {
            return false;
        }
    }
    return true;
}

math::linal::SparseMatrix::value_type math::linal::SparseMatrix::get(size_t row, size_t col) const {
    if (row >= get_height() || col >= get_width()) {
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    }
    std::pair<size_t, size_t> key = { row, col };
    if (auto it = m_data.find(key); it != m_data.end()) {
        return it->second;
    }
    return 0;
}

void math::linal::SparseMatrix::set(size_t row, size_t col, value_type value) {
    if (row >= get_height() || col >= get_width()) {
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    }
    if (dcmp(value) != 0) {
        m_data[{row, col}] = value;
    }    
}

size_t math::linal::SparseMatrix::get_width() const {
    return m_width;
}

size_t math::linal::SparseMatrix::get_height() const {
    return m_height;
}

void math::linal::SparseMatrix::shrink() {
    for (auto it = m_data.begin(); it != m_data.end();) {
        if (dcmp(it->second) == 0) {
            it = m_data.erase(it);
        }
        else {
            ++it;
        }
    }
}

math::linal::SparseMatrix::value_type math::linal::SparseMatrix::det() const {
    value_type determinant = 1.0;
    if (is_empty())
        return determinant;

    // ...
    return determinant;
}

std::vector<std::pair<size_t, math::linal::SparseMatrix::complex_value_type>> math::linal::SparseMatrix::get_eigenvalues() const {
    std::vector<std::pair<size_t, complex_value_type>> eigenvalues;
    // ...
    return eigenvalues;
}

std::vector<std::pair<math::linal::SparseMatrix::complex_value_type, std::vector<math::linal::DVector>>> math::linal::SparseMatrix::get_eigenvectors() const {
    std::vector<std::pair<complex_value_type, std::vector<math::linal::DVector>>> eigenvectors;
    auto eigenvalues = get_eigenvalues();
    // ...
    return eigenvectors;
}

math::linal::SparseMatrix math::linal::SparseMatrix::identity_matrix(size_t n) {
    SparseMatrix res(n, n);
    for (size_t i = 0; i < n; ++i) {
        res.m_data.insert({ { i, i }, 1.0 });
    }
    return res;
}

math::linal::SparseMatrix math::linal::SparseMatrix::elementary_matrix_unit(size_t n, size_t m, size_t i, size_t j, double value) {
    SparseMatrix res(n, m);
    res.set(i, j, value);
    return res;
}

void math::linal::SparseMatrix::swap(SparseMatrix& other) noexcept {
    std::swap(m_width, other.m_width);
    std::swap(m_height, other.m_height);
    std::swap(m_data, other.m_data);
}

void math::linal::swap(SparseMatrix& m1, SparseMatrix& m2) noexcept {
    m1.swap(m2);
}

bool math::linal::operator==(const SparseMatrix& m1, const SparseMatrix& m2) {
    if (&m1 == &m2)
        return true;
    if (m1.get_height() != m2.get_height() || m1.get_width() != m2.get_width())
        return false;

    for (const auto& [key, value] : m1.m_data) {
        auto it = m2.m_data.find(key);
        if (it == m2.m_data.end()) {
            if (dcmp(value) != 0) {
                return false;
            }
        }
        else if (dcmp(value, it->second) != 0) {
            return false;
        }
    }

    for (const auto& [key, value] : m2.m_data) {
        if (auto it = m1.m_data.find(key); it == m1.m_data.end() && dcmp(value) != 0) {
            return false;
        }
    }
    return true;
}

bool math::linal::operator!=(const SparseMatrix& m1, const SparseMatrix& m2) {
    return !(m1 == m2);
}

math::linal::SparseMatrix math::linal::operator*(const SparseMatrix& matrix, SparseMatrix::value_type scalar) {
    SparseMatrix res(matrix);
    res *= scalar;
    return res;
}

math::linal::SparseMatrix math::linal::operator*(SparseMatrix::value_type scalar, const SparseMatrix& matrix) {
    return matrix * scalar;
}

math::linal::SparseMatrix math::linal::operator/(const SparseMatrix& matrix, SparseMatrix::value_type scalar) {
    return matrix * (1. / scalar);
}

math::linal::DVector math::linal::operator*(const SparseMatrix& matrix, const DVector& vector) {
    if (matrix.get_width() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }

    DVector res(matrix.get_height(), 0.0);
    for (const auto& [key, value] : matrix.m_data) {
        res[key.first] += value * vector[key.second];
    }
    return res;
}

math::linal::DVector math::linal::operator*(const DVector& vector, const SparseMatrix& matrix) {
    if (matrix.get_height() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }

    DVector res(matrix.get_width(), 0.0);
    for (const auto& [key, value] : matrix.m_data) {
        res[key.second] += value * vector[key.first];
    }
    return res;
}

math::linal::SparseMatrix math::linal::operator*(const SparseMatrix& m1, const SparseMatrix& m2) {
    if (m1.get_width() != m2.get_height()) {
        throw std::invalid_argument("The matrices must be the same size");
    }

    const size_t n = m1.get_height();  
    const size_t p = m2.get_width();   

    SparseMatrix res(n, p);

    // Для эффективного доступа к столбцам m2 предварительно группируем их
    std::vector<std::unordered_map<size_t, SparseMatrix::value_type>> m2_cols(p);
    for (const auto& [key, val] : m2.m_data) {
        m2_cols[key.second][key.first] = val;
    }

    for (const auto& [key1, val1] : m1.m_data) {
        const auto& [i, k] = key1;  
        for (const auto& [j, val2] : m2_cols[k]) {
            res.m_data[{i, j}] += val1 * val2; 
        }
    }

    return res;
}

math::linal::SparseMatrix math::linal::operator+(const SparseMatrix& m1, const SparseMatrix& m2) {
    if (m1.get_height() != m2.get_height() || m1.get_width() != m2.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }
    SparseMatrix res(m1);
    res += m2;
    return res;
}

math::linal::SparseMatrix math::linal::operator-(const SparseMatrix& m1, const SparseMatrix& m2) {
    if (m1.get_height() != m2.get_height() || m1.get_width() != m2.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }
    SparseMatrix res(m1);
    res -= m2;
    return res;
}

math::linal::SparseMatrix math::linal::operator-(const SparseMatrix& matrix) {
    return -1 * matrix;
}

math::linal::SparseMatrix math::linal::inversed(const SparseMatrix& matrix) {
    if (matrix.get_width() == 0) {
        throw std::invalid_argument("Matrix is empty");
    }
    
    SparseMatrix res(matrix.get_width(), matrix.get_height());
    // ...
    return res;
}

math::linal::SparseMatrix math::linal::transposed(const SparseMatrix& matrix) {
    SparseMatrix res(matrix.get_height(), matrix.get_width());

    for (const auto& [key, val] : matrix.m_data) {
        res.m_data.insert({ {key.second, key.first}, val });
    }
    return res;
}
