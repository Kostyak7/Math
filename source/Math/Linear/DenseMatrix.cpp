#include <Utils/Math/Linear/DenseMatrix.h>

#include <Utils/Math/Common.h>

#include <stdexcept>

math::linal::DenseMatrix::DenseMatrix() noexcept
    : IMatrixFromVector()
    , m_height(0)
    , m_width(0)
{
}

math::linal::DenseMatrix::DenseMatrix(size_t height, size_t width, const value_type& default_value)
    : IMatrixFromVector(height * width, default_value)
    , m_height(height)
    , m_width(width)
{
}

math::linal::DenseMatrix::DenseMatrix(const DenseMatrix& matrix) 
    : IMatrixFromVector(matrix.m_data)
    , m_height(matrix.m_height)
    , m_width(matrix.m_width)

{
}

math::linal::DenseMatrix::DenseMatrix(DenseMatrix&& matrix) noexcept
    : IMatrixFromVector(std::move(matrix.m_data))
    , m_height(matrix.m_height)
    , m_width(matrix.m_width)
{
}

math::linal::DenseMatrix::DenseMatrix(std::initializer_list<value_type> list) 
    : IMatrixFromVector(list.size())
    , m_height(list.size())
    , m_width(1)
{
    auto it = list.begin();
    for (size_t i = 0; i < m_height; ++i, ++it) {
        _(i, 0) = *it;
    }
}

math::linal::DenseMatrix::DenseMatrix(std::initializer_list<std::initializer_list<value_type>> list) 
    : IMatrixFromVector(list.size() * list.begin()->size())
    , m_height(list.size())
    , m_width(list.begin()->size())
{
    size_t width = 0;
    for (const auto& lst : list) {
        width = std::max(width, lst.size());
    }
    if (width != m_width) {
        reshape(m_height, width);
    }

    auto it = list.begin();
    for (size_t r = 0; r < m_height; ++r, ++it) {
        size_t c = 0;
        for (const auto& el : *it) {
            _(r, c) = el;
            ++c;
            m_data;
        }
    }
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator=(const DenseMatrix& matrix) {
    if (this == &matrix) return *this;
    m_height = matrix.m_height;
    m_width = matrix.m_width;
    m_data = matrix.m_data;
    return *this;
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator=(DenseMatrix&& matrix) noexcept {
    if (this == &matrix) return *this;
    m_height = matrix.m_height;
    m_width = matrix.m_width;
    m_data = std::move(matrix.m_data);
    return *this;
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator*=(value_type scalar) {
    for (auto& val : m_data) {
        val *= scalar;
    }
    return *this;
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator/=(value_type scalar) {
    return *this *= (1. / scalar);
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator+=(const DenseMatrix& other) {
    if (get_width() != other.get_width() || get_height() != other.get_height())
        throw std::invalid_argument("The matrices must be the same size");
    for (size_t i = 0; i < m_data.size(); ++i) {
        m_data[i] += other.m_data[i];
    }
    return *this;
}

math::linal::DenseMatrix& math::linal::DenseMatrix::operator-=(const DenseMatrix& other) {
    if (get_width() != other.get_width() || get_height() != other.get_height())
        throw std::invalid_argument("The matrices must be the same size");
    for (size_t i = 0; i < m_data.size(); ++i) {
        m_data[i] -= other.m_data[i];
    }
    return *this;
}

void math::linal::DenseMatrix::swap(DenseMatrix& matrix) noexcept {
    std::swap(m_height, matrix.m_height);
    std::swap(m_width, matrix.m_width);
    std::swap(m_data, matrix.m_data);
}

void math::linal::DenseMatrix::reshape(size_t height, size_t width) {
    m_height = height;
    m_width = width;
    m_data.resize(m_height * m_width);
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
                sum += L._(i, k) * L._(j, k);
            }

            if (i == j) {
                value_type diag = _(i, j) - sum;
                if (diag <= 0) 
                    return false; 

                L._(i, i) = std::sqrt(diag);
            }
            else {
                L._(i, j) = (_(i, j) - sum) / L._(j, j);
            }
        }
    }
    return true;
}

math::linal::DenseMatrix::value_type math::linal::DenseMatrix::get(size_t row, size_t col) const {
    if (row >= get_height() || col >= get_width())
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    return _(row, col);
}

void math::linal::DenseMatrix::set(size_t row, size_t col, value_type value) {
    if (row >= get_height() || col >= get_width())
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    _(row, col) = value;
}

size_t math::linal::DenseMatrix::get_width() const {
    return m_width;
}

size_t math::linal::DenseMatrix::get_height() const {
    return m_height;
}

math::linal::DenseMatrix::value_type math::linal::DenseMatrix::det() const {
    if (!is_square()) {
        throw std::runtime_error("Determinant is defined only for square matrices");
    }

    size_t n = get_width();
    if (n == 0) return 0.0;
    if (n == 1) return _(0, 0);
    if (n == 2) return _(0, 0) * _(1, 1) - _(0, 1) * _(1, 0);

    DenseMatrix temp(*this);
    value_type det = 1.0;

    for (size_t i = 0; i < n; ++i) {
        size_t pivot = i;
        for (size_t j = i + 1; j < n; ++j) {
            if (fabs(temp._(j, i)) > fabs(temp._(pivot, i))) {
                pivot = j;
            }
        }

        if (pivot != i) {
            temp.swap_rows(i, pivot);
            det = -det;
        }

        if (temp._(i, i) == 0) 
            return 0.0;

        det *= temp._(i, i);

        for (size_t j = i + 1; j < n; ++j) {
            value_type factor = temp._(j, i) / temp._(i, i);
            for (size_t k = i + 1; k < n; ++k) {
                temp._(j, k) -= factor * temp._(i, k);
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

void math::linal::DenseMatrix::swap_rows(size_t r1, size_t r2) {
    value_type tmp;
    auto row1 = get_row(r1);
    for (auto it1 = row1.begin(), it2 = get_row(r2).begin(); it1 != row1.end(); ++it1, ++it2) {
        tmp = *it1;
        *it1 = *it2;
        *it2 = tmp;
    }
}

void math::linal::DenseMatrix::swap_columns(size_t c1, size_t c2) {
    value_type tmp;
    for (size_t j = 0; j < m_height; ++j) {
        tmp = _(j, c1);
        _(j, c1) = _(j, c2);
        _(j, c2) = tmp;
    }
}

math::linal::DenseMatrix math::linal::DenseMatrix::identity_matrix(size_t n) {
    DenseMatrix res(n, n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        res._(i, i) = 1.0;
    }
    return res;
}

math::linal::DenseMatrix math::linal::DenseMatrix::elementary_matrix_unit(size_t n, size_t m, size_t i, size_t j, double value) {
    DenseMatrix res(n, m);
    res.set(i, j, value);
    return res;
}

math::linal::DenseMatrix::ConstProxyVector math::linal::DenseMatrix::get_row(size_t row) const {
    return ConstProxyVector(m_width, std::next(m_data.begin(), row * m_width));
}

math::linal::DenseMatrix::ProxyVector math::linal::DenseMatrix::get_row(size_t row) {
    return ProxyVector(m_width, std::next(m_data.begin(), row * m_width));
}

bool math::linal::DenseMatrix::check_row_index(size_t row) const {
    return row < m_height;
}

size_t math::linal::DenseMatrix::calculate_index(size_t row, size_t col) const {
    return col + row * m_width;
}

void math::linal::swap(DenseMatrix& m1, DenseMatrix& m2) noexcept {
    m1.swap(m2);
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

math::linal::DenseMatrix math::linal::operator*(DenseMatrix matrix, DenseMatrix::value_type scalar) {
    matrix *= scalar;
    return matrix;
}

math::linal::DenseMatrix math::linal::operator*(DenseMatrix::value_type scalar, DenseMatrix matrix) {
    return matrix * scalar;
}

math::linal::DenseMatrix math::linal::operator/(DenseMatrix matrix, DenseMatrix::value_type scalar) {
    return matrix * (1. / scalar);
}

math::linal::DVector math::linal::operator*(const DenseMatrix& matrix, const DVector& vector) {
    if (matrix.get_width() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }
    const size_t n = matrix.get_height();
    const size_t m = matrix.get_width();
    DVector res(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        const auto row = matrix[i];
        for (size_t j = 0; j < m; ++j) {
            res[i] = row[j] * vector[j];
        }
    }
    return res;
}

math::linal::DVector math::linal::operator*(const DVector& vector, const DenseMatrix& matrix) {
    if (matrix.get_height() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }
    const size_t n = matrix.get_height();
    const size_t m = matrix.get_width();
    DVector res(m, 0.0);
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            res[i] += vector[j] * matrix[j][i];
        }
    }
    return res;
}

math::linal::DenseMatrix math::linal::operator*(const DenseMatrix& m1, const DenseMatrix& m2) {
    if (m1.get_width() != m2.get_height()) {
        throw std::invalid_argument("The matrices must be the same size");
    }
    DenseMatrix res(m1.get_height(), m2.get_width());
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
    if (!matrix.is_square()) {
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
