#include <Utils/Math/Linear/BandMatrix.h>

#include <Utils/Math/Common.h>

#include <stdexcept>

math::linal::BandMatrix::BandMatrix() noexcept
    : IMatrixFromVector()
    , m_size(0)
    , m_isl(0)
{
}

math::linal::BandMatrix::BandMatrix(size_t size, size_t isl, const value_type& default_value)
    : IMatrixFromVector(size* isl, default_value)
    , m_size(size)
    , m_isl(isl)
{
}

math::linal::BandMatrix::BandMatrix(const BandMatrix& matrix)
    : IMatrixFromVector(matrix.m_data)
    , m_size(matrix.m_size)
    , m_isl(matrix.m_isl)

{
}

math::linal::BandMatrix::BandMatrix(BandMatrix&& matrix) noexcept
    : IMatrixFromVector(std::move(matrix.m_data))
    , m_size(matrix.m_size)
    , m_isl(matrix.m_isl)
{
}

math::linal::BandMatrix::BandMatrix(std::initializer_list<value_type> list)
    : IMatrixFromVector(list.size())
    , m_size(list.size())
    , m_isl(1)
{
    auto it = list.begin();
    for (size_t i = 0; i < m_size; ++i, ++it) {
        _(i, 0) = *it;
    }
}

math::linal::BandMatrix::BandMatrix(std::initializer_list<std::initializer_list<value_type>> list)
    : IMatrixFromVector(list.size() * list.begin()->size())
    , m_size(list.size())
    , m_isl(list.begin()->size())
{
    size_t isl = 0;
    for (const auto& lst : list) {
        isl = std::max(isl, lst.size());
    }
    if (isl != m_isl) {
        reshape(m_size, isl);
    }

    auto rows_it = list.begin();
    for (size_t r = 0; r < m_size; ++r, ++rows_it) {
        auto it = rows_it->begin();
        for (size_t c = 0; c < std::min(m_isl, rows_it->size()); ++c, ++it) {
            _(r, c) = *it;
        }        
    }
}

math::linal::BandMatrix& math::linal::BandMatrix::operator=(const BandMatrix& matrix) {
    if (this == &matrix) return *this;
}

math::linal::BandMatrix& math::linal::BandMatrix::operator=(BandMatrix&& matrix) noexcept {
    if (this == &matrix) return *this;
    m_size = matrix.m_size;
    m_isl = matrix.m_isl;
    m_data = std::move(matrix.m_data);
    return *this;
}

math::linal::BandMatrix& math::linal::BandMatrix::operator*=(value_type scalar) {
    for (auto& val : m_data) {
        val *= scalar;
    }
    return *this;
}

math::linal::BandMatrix& math::linal::BandMatrix::operator/=(value_type scalar) {
    return *this *= (1. / scalar);
}

math::linal::BandMatrix& math::linal::BandMatrix::operator+=(const BandMatrix& other) {
    if (get_height() != other.get_height() || get_width() != other.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }

    if (m_isl < other.m_isl) {
        set(0, other.m_isl - 1, 0.);
    }
    for (size_t r = 0, isl = std::min(m_isl, other.m_isl); r < m_size; ++r) {
        for (size_t c = 0; c < isl; ++c) {
            _(r, c) += other._(r, c);
        }
    }
    return *this;
}

math::linal::BandMatrix& math::linal::BandMatrix::operator-=(const BandMatrix& other) {
    if (get_height() != other.get_height() || get_width() != other.get_width()) {
        throw std::invalid_argument("The matrices must be the same size");
    }

    if (m_isl < other.m_isl) {
        set(0, other.m_isl - 1, 0.);
    }
    for (size_t r = 0, isl = std::min(m_isl, other.m_isl); r < m_size; ++r) {
        for (size_t c = 0; c < isl; ++c) {
            _(r, c) -= other._(r, c);
        }
    }
    return *this;
}

void math::linal::BandMatrix::swap(BandMatrix& matrix) noexcept {
    std::swap(m_size, matrix.m_size);
    std::swap(m_isl, matrix.m_isl);
    std::swap(m_data, matrix.m_data);
}

void math::linal::BandMatrix::reshape(size_t size, size_t isl) {
    m_size = size;
    m_isl = isl;
    m_data.resize(m_size * m_isl);
}

bool math::linal::BandMatrix::is_zero() const {
    for (const auto& val : m_data) {
        if (dcmp(val) != 0) {
            return false;
        }
    }
    return true;
}

bool math::linal::BandMatrix::is_identity() const {
    if (!is_square())
        return false;
    for (size_t i = 0; i < m_data.size();) {
        if (dcmp(m_data[i], 1.) != 0) return false;
        ++i;
        for (size_t j = 1; j < m_isl; ++i, ++j) {
            if (dcmp(m_data[i]) != 0)
                return false;
        }
    }
    return true;
}

bool math::linal::BandMatrix::is_diagonal() const {
    if (!is_square())
        return false;
    for (size_t i = 0; i < m_data.size();) {
        if (dcmp(m_data[i]) == 0) return false;
        ++i;
        for (size_t j = 1; j < m_isl; ++i, ++j) {
            if (dcmp(m_data[i]) != 0)
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
    if (col >= m_isl + row || row >= m_isl + col) {
        return 0.0;
    }
    if (col < row) {
        return _(col, row - col);
    }
    return _(row, col - row);
}

void math::linal::BandMatrix::set(size_t row, size_t col, value_type value) {
    if (row >= get_height() || col >= get_width()) {
        throw std::invalid_argument("Going beyond the boundaries of the matrix");
    }
    if (col >= m_isl + row) {
        set_isl(col - row + 1);
    }
    else if (row >= m_isl + col) {
        set_isl(row - col + 1);
    }
    if (col < row) {
        _(col, row - col) = value;
    }
    _(row, col - row) = value;
}

size_t math::linal::BandMatrix::get_width() const {
    return m_size;
}

size_t math::linal::BandMatrix::get_height() const {
    return m_size;
}

size_t math::linal::BandMatrix::get_isl() const {
    return m_isl;
}

void math::linal::BandMatrix::set_isl(size_t width) {
    if (width == m_isl) return;
    std::vector<value_type> new_data(m_size * width, 0.0);
    for (size_t i = 0, j = 0; j < m_size; ++j) {
        size_t offset = m_isl * j;
        size_t new_offset = width * j;
        for (; i < std::min(width, m_isl); ++i) {
            new_data[i + new_offset] = m_data[i + offset];
        }
        for (; i < width; ++i) {
            new_data[i + new_offset] = 0.0;
        }
    }    
    m_isl = width;
    m_data = std::move(new_data);
}

void math::linal::BandMatrix::shrink() {
    if (m_isl == 0)
        return;
    const size_t old_isl = m_isl;
    size_t new_isl = 1;
    for (size_t row = 0; row < get_height(); ++row) {
        size_t zero_count = 0;
        for (size_t col = old_isl - 1; col >= 0; --col) {
            if (dcmp(_(row, col)) != 0) {
                break;
            }
            ++zero_count;
        }
        if (zero_count == 0)
            return;
        new_isl = std::max(new_isl, old_isl - zero_count);
    }
    set_isl(new_isl);
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

std::vector<std::pair<math::linal::BandMatrix::complex_value_type, std::vector<math::linal::DVector>>> math::linal::BandMatrix::get_eigenvectors() const {
    std::vector<std::pair<complex_value_type, std::vector<math::linal::DVector>>> eigenvectors;
    auto eigenvalues = get_eigenvalues();
    // ...
    return eigenvectors;
}

math::linal::BandMatrix math::linal::BandMatrix::identity_matrix(size_t n, size_t isl) {
    BandMatrix res(n, isl, 0.0);
    for (size_t i = 0; i < isl; ++i) {
        res._(i, 0) = 1.0;
    }
    return res;
}

math::linal::BandMatrix math::linal::BandMatrix::elementary_matrix_unit(size_t n, size_t i, size_t j, value_type value) {
    BandMatrix res(n, i);
    res.set(i, j, value);
    return res;
}

math::linal::BandMatrix::ConstProxyVector math::linal::BandMatrix::get_row(size_t row) const {
    return { m_isl, std::next(m_data.begin(), row * m_isl) };
}

math::linal::BandMatrix::ProxyVector math::linal::BandMatrix::get_row(size_t row) {
    return { m_isl, std::next(m_data.begin(), row * m_isl) };
}

bool math::linal::BandMatrix::check_row_index(size_t row) const {
    return row < m_size;
}

size_t math::linal::BandMatrix::calculate_index(size_t row, size_t col) const {
    return col + row * m_isl;
}

void math::linal::swap(BandMatrix& m1, BandMatrix& m2) noexcept {
    m1.swap(m2);
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
            const DVector& vec = (m1.get_isl() < m2.get_isl()) ? m2[row] : m1[row];
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
    return matrix * scalar;
}

math::linal::BandMatrix math::linal::operator/(const BandMatrix& matrix, BandMatrix::value_type scalar) {
    return matrix * (1. / scalar);
}

math::linal::DVector math::linal::operator*(const BandMatrix& matrix, const DVector& vector) {
    if (matrix.get_width() != vector.size()) {
        throw std::invalid_argument("The matrix and the vector must be the same size");
    }

    size_t n = matrix.get_height();
    size_t w = matrix.get_isl();
    DVector res(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = std::max(0, static_cast<int>(i) - static_cast<int>(w)); j <= std::min(n - 1, i + w); ++j) {
            res[i] += matrix.get(i, j) * vector[j];
        }
    }
    return res;
}

math::linal::DVector math::linal::operator*(const DVector& vector, const BandMatrix& matrix) {
    return matrix * vector;
}

math::linal::BandMatrix math::linal::operator*(const BandMatrix& m1, const BandMatrix& m2) {
    if (m1.get_width() != m2.get_height()) {
        throw std::invalid_argument("The matrices must be the same size");
    }

    size_t n = m1.get_height();
    size_t w = m1.get_isl();  
    BandMatrix result(n, w + 1, 0.0); 

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
    return -1. * matrix;
}

math::linal::BandMatrix math::linal::inversed(const BandMatrix& matrix) {
    if (matrix.is_empty()) {
        throw std::invalid_argument("Matrix is empty");
    }

    size_t n = matrix.get_height();
    BandMatrix L(n, n, 0.0);
    // ...
    BandMatrix res(n, n, 0.0);
    // ...
    return res;
}
