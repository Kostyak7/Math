#include <Utils/Math/Linear/DVector.h>

#include <Utils/Math/Common.h>

#include <stdexcept>

math::linal::DVector::iterator math::linal::DVector::begin() noexcept {
    return m_data.begin();
}

math::linal::DVector::const_iterator math::linal::DVector::begin() const noexcept {
    return m_data.begin();
}

math::linal::DVector::iterator math::linal::DVector::end() noexcept {
    return m_data.end();
}

math::linal::DVector::const_iterator math::linal::DVector::end() const noexcept {
    return m_data.end();
}

math::linal::DVector::reverse_iterator math::linal::DVector::rbegin() noexcept {
    return m_data.rbegin();
}

math::linal::DVector::const_reverse_iterator math::linal::DVector::rbegin() const noexcept {
    return m_data.rbegin();
}

math::linal::DVector::reverse_iterator math::linal::DVector::rend() noexcept {
    return m_data.rend();
}

math::linal::DVector::const_reverse_iterator math::linal::DVector::rend() const noexcept {
    return m_data.rend();
}

math::linal::DVector::const_iterator math::linal::DVector::cbegin() const noexcept {
    return m_data.cbegin();
}

math::linal::DVector::const_iterator math::linal::DVector::cend() const noexcept {
    return m_data.cend();
}

math::linal::DVector::const_reverse_iterator math::linal::DVector::crbegin() const noexcept {
    return m_data.crbegin();
}

math::linal::DVector::const_reverse_iterator math::linal::DVector::crend() const noexcept {
    return m_data.crend();
}

const math::linal::DVector::value_type& math::linal::DVector::front() const noexcept {
    return m_data.front();
}

math::linal::DVector::value_type& math::linal::DVector::front() noexcept {
    return m_data.front();
}

const math::linal::DVector::value_type& math::linal::DVector::back() const noexcept {
    return m_data.back();
}

math::linal::DVector::value_type& math::linal::DVector::back() noexcept {
    return m_data.back();
}

math::linal::DVector::DVector(size_t size, value_type default_val)
    : m_data(size, default_val)
{
}

math::linal::DVector::DVector(std::initializer_list<value_type> init)
    : m_data(init)
{
}

math::linal::DVector::DVector(const base_container_type& container)
    : m_data(container)
{
}

math::linal::DVector::DVector(const DVector& other)
    : m_data(other.m_data)
{
}

math::linal::DVector::DVector(DVector&& other) noexcept
    : m_data(std::move(other.m_data))
{
}

math::linal::DVector& math::linal::DVector::operator=(const DVector& other) {
    if (this == &other) return *this;
    m_data = other.m_data;
    return *this;
}

math::linal::DVector& math::linal::DVector::operator=(DVector&& other) noexcept {
    if (this == &other) return *this;
    m_data = std::move(other.m_data);
    return *this;
}

math::linal::DVector& math::linal::DVector::operator*=(value_type scalar) {
    for (auto& val : *this) {
        val *= scalar;
    }
    return *this;
}

math::linal::DVector& math::linal::DVector::operator/=(value_type scalar) {
    if (dcmp(scalar) == 0) {
        throw std::invalid_argument("Impossible to divide by zero");
    }
    return *this *= (1.0 / scalar);
}

math::linal::DVector& math::linal::DVector::operator+=(const DVector& other) {
    if (size() != other.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] += other[i];
    }
    return *this;
}

math::linal::DVector& math::linal::DVector::operator-=(const DVector& other) {
    if (size() != other.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] -= other[i];
    }
    return *this;
}

void math::linal::DVector::swap(DVector& other) noexcept {
    std::swap(m_data, other.m_data);
}

size_t math::linal::DVector::size() const {
    return m_data.size();
}

bool math::linal::DVector::empty() const {
    return m_data.empty();
}

void math::linal::DVector::resize(size_t new_size) {
    m_data.resize(new_size);
}

void math::linal::DVector::resize(size_t new_size, value_type val) {
    m_data.resize(new_size, val);
}

void math::linal::DVector::fill(value_type value) {
    for (auto& el : m_data) {
        el = value;
    }
}

bool math::linal::DVector::is_zero() const {
    for (const auto& val : *this) {
        if (dcmp(val) != 0) {
            return false;
        }
    }
    return true;
}

const math::linal::DVector::value_type& math::linal::DVector::operator[](size_t pos) const noexcept {
    return m_data[pos];
}

math::linal::DVector::value_type& math::linal::DVector::operator[](size_t pos) noexcept {
    return m_data[pos];
}

const math::linal::DVector::value_type& math::linal::DVector::at(size_t pos) const {
    return m_data.at(pos);
}

math::linal::DVector::value_type& math::linal::DVector::at(size_t pos) {
    return m_data.at(pos);
}

math::linal::DVector::value_type math::linal::DVector::dot(const DVector& other) const {
    if (size() != other.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    value_type result = 0.0;
    for (size_t i = 0; i < size(); ++i) {
        result += (*this)[i] * other[i];
    }
    return result;
}

math::linal::DVector::value_type math::linal::DVector::norm() const {
    return std::sqrt(dot(*this));
}

math::linal::DVector::value_type math::linal::DVector::length() const {
    return norm();
}

math::linal::DVector& math::linal::DVector::normalize() {
    value_type norm = this->norm();
    if (dcmp(norm) == 0) {
        throw std::runtime_error("Cannot normalize zero vector");
    }
    return *this /= norm;
}

void math::linal::swap(DVector& v1, DVector& v2) noexcept {
    v1.swap(v2);
}

bool math::linal::operator==(const DVector& v1, const DVector& v2) {
    if (v1.size() != v2.size()) 
        return false;
    for (size_t i = 0; i < v1.size(); ++i) {
        if (dcmp(v1[i], v2[i]) != 0) 
            return false;
    }
    return true;
}

bool math::linal::operator!=(const DVector& v1, const DVector& v2) {
    return !(v1 == v2);
}

math::linal::DVector math::linal::operator*(DVector vector, DVector::value_type scalar) {
    vector *= scalar;
    return vector;
}

math::linal::DVector math::linal::operator*(DVector::value_type scalar, DVector vector) {
    return vector * scalar;
}

math::linal::DVector math::linal::operator/(DVector vector, DVector::value_type scalar) {
    return vector * (1. / scalar);
}

math::linal::DVector math::linal::operator+(const DVector& v1, const DVector& v2) {
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    DVector result(v1);
    return result += v2;
}

math::linal::DVector math::linal::operator-(const DVector& v1, const DVector& v2) {
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    DVector result(v1);
    return result -= v2;
}

math::linal::DVector math::linal::operator-(const DVector& vector) {
    return -1 * vector;
}

math::linal::DVector math::linal::normalized(DVector vector) {
    DVector::value_type norm = 0.0;
    for (const auto& val : vector) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    if (dcmp(norm) == 0) {
        throw std::runtime_error("Cannot normalize zero vector");
    }
    return vector * (1.0 / norm);
}
