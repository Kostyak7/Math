#include <Utils/Math/Linear/IMatrixFromVector.h>

#include <Utils/Math/Common.h>

#include <stdexcept>

//////////////////////////////////////////
//          IMatrixFromVector           //
//////////////////////////////////////////

math::linal::IMatrixFromVector::IMatrixFromVector(size_t size, const value_type& default_value) 
    : m_data(size, default_value)
{
}

math::linal::IMatrixFromVector::IMatrixFromVector(std::vector<value_type> data) 
    : m_data(std::move(data))
{
}

math::linal::IMatrixFromVector::ConstProxyVector math::linal::IMatrixFromVector::operator[](size_t row) const {
    return get_row(row);
}

math::linal::IMatrixFromVector::ProxyVector math::linal::IMatrixFromVector::operator[](size_t row) {
    return get_row(row);
}

math::linal::IMatrixFromVector::ConstProxyVector math::linal::IMatrixFromVector::at(size_t row) const {
    if (!check_row_index(row)) {
        throw std::invalid_argument("Row index out of range!");
    }
    return get_row(row);
}

math::linal::IMatrixFromVector::ProxyVector math::linal::IMatrixFromVector::at(size_t row) {
    if (!check_row_index(row)) {
        throw std::invalid_argument("Row index out of range!");
    }
    return get_row(row);
}

void math::linal::IMatrixFromVector::clear() {
    m_data.clear();
}

const math::linal::IMatrixFromVector::value_type& math::linal::IMatrixFromVector::_(size_t row, size_t col) const {
    return m_data[calculate_index(row, col)];
}

math::linal::IMatrixFromVector::value_type& math::linal::IMatrixFromVector::_(size_t row, size_t col) {
    return m_data[calculate_index(row, col)];
}

//////////////////////////////////////////
//           ConstProxyVector           //
//////////////////////////////////////////

math::linal::IMatrixFromVector::ConstProxyVector::const_iterator math::linal::IMatrixFromVector::ConstProxyVector::begin() const {
    return m_front;
}

math::linal::IMatrixFromVector::ConstProxyVector::const_iterator math::linal::IMatrixFromVector::ConstProxyVector::end() const {
    return std::next(m_front, m_size);
}

math::linal::IMatrixFromVector::ConstProxyVector::ConstProxyVector(size_t size, const_iterator front)
    : m_size(size)
    , m_front(front)
{
}

math::linal::IMatrixFromVector::ConstProxyVector::operator math::linal::DVector() const {
    DVector res(m_size);
    auto it = m_front;
    for (size_t i = 0; i < m_size; ++i, ++it) {
        res[i] = *it;
    }
    return res;
}

math::linal::IMatrixFromVector::value_type math::linal::IMatrixFromVector::ConstProxyVector::operator[](size_t col) const {
    return *std::next(m_front, col);
}

math::linal::IMatrixFromVector::value_type math::linal::IMatrixFromVector::ConstProxyVector::at(size_t col) const {
    if (col >= m_size) {
        throw std::invalid_argument("ProxyVector out of range!");
    }
    return *std::next(m_front, col);
}

size_t math::linal::IMatrixFromVector::ConstProxyVector::size() const {
    return m_size;
}

bool math::linal::IMatrixFromVector::ConstProxyVector::empty() const {
    return m_size == 0;
}

bool math::linal::operator==(const IMatrixFromVector::ConstProxyVector& v1, const IMatrixFromVector::ConstProxyVector& v2) {
    if (&v1 == &v2)
        return true;
    if (v1.size() != v2.size())
        return false;
    for (auto it1 = v1.begin(), it2 = v2.begin(); it1 != v1.end(); ++it1, ++it2) {
        if (*it1 != *it2) {
            return false;
        }
    }
    return true;
}

bool math::linal::operator!=(const IMatrixFromVector::ConstProxyVector& v1, const IMatrixFromVector::ConstProxyVector& v2) {
    return !(v1 == v2);
}

//////////////////////////////////////////
//             ProxyVector              //
//////////////////////////////////////////

math::linal::IMatrixFromVector::ProxyVector::const_iterator math::linal::IMatrixFromVector::ProxyVector::begin() const {
    return m_front;
}

math::linal::IMatrixFromVector::ProxyVector::const_iterator math::linal::IMatrixFromVector::ProxyVector::end() const {
    return std::next(m_front, m_size); 
}

math::linal::IMatrixFromVector::ProxyVector::iterator math::linal::IMatrixFromVector::ProxyVector::begin() {
    return m_front;
}

math::linal::IMatrixFromVector::ProxyVector::iterator math::linal::IMatrixFromVector::ProxyVector::end() {
    return std::next(m_front, m_size);
}

math::linal::IMatrixFromVector::ProxyVector::ProxyVector(size_t size, iterator front)
    : m_size(size)
    , m_front(front)
{
}

math::linal::IMatrixFromVector::ProxyVector::operator math::linal::DVector() const {
    DVector res(m_size);
    auto it = m_front;
    for (size_t i = 0; i < m_size; ++i, ++it) {
        res[i] = *it;
    }
    return res;
}

math::linal::IMatrixFromVector::value_type math::linal::IMatrixFromVector::ProxyVector::operator[](size_t col) const {
    return *std::next(m_front, col);
}

math::linal::IMatrixFromVector::value_type& math::linal::IMatrixFromVector::ProxyVector::operator[](size_t col) {
    return *std::next(m_front, col);
}

math::linal::IMatrixFromVector::value_type math::linal::IMatrixFromVector::ProxyVector::at(size_t col) const {
    if (col >= m_size) {
        throw std::invalid_argument("ProxyVector out of range!");
    }
    return *std::next(m_front, col);
}

math::linal::IMatrixFromVector::value_type& math::linal::IMatrixFromVector::ProxyVector::at(size_t col) {
    if (col >= m_size) {
        throw std::invalid_argument("ProxyVector out of range!");
    }
    return *std::next(m_front, col);
}

size_t math::linal::IMatrixFromVector::ProxyVector::size() const {
    return m_size;
}

bool math::linal::IMatrixFromVector::ProxyVector::empty() const {
    return m_size == 0;
}

bool math::linal::operator==(const IMatrixFromVector::ProxyVector& v1, const IMatrixFromVector::ProxyVector& v2) {
    if (&v1 == &v2) 
        return true;
    if (v1.size() != v2.size()) 
        return false;
    for (auto it1 = v1.begin(), it2 = v2.begin(); it1 != v1.end(); ++it1, ++it2) {
        if (*it1 != *it2) {
            return false;
        }
    }
    return true;
}

bool math::linal::operator!=(const IMatrixFromVector::ProxyVector& v1, const IMatrixFromVector::ProxyVector& v2) {
    return !(v1 == v2);
}