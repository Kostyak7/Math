#include "FVector.h"

#include <Utils/Math/Common.h>

#include <stdexcept>

math::linal::FVector& math::linal::FVector::operator*=(value_type scalar) {
    for (auto& val : *this) {
        val *= scalar;
    }
    return *this;
}

math::linal::FVector& math::linal::FVector::operator/=(value_type scalar) {
    if (dcmp(scalar) == 0) {
        throw std::invalid_argument("Impossible to divide by zero");
    }
    return *this *= (1.0 / scalar);
}

math::linal::FVector& math::linal::FVector::operator+=(const FVector& other) {
    if (size() != other.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] += other[i];
    }
    return *this;
}

math::linal::FVector& math::linal::FVector::operator-=(const FVector& other) {
    if (size() != other.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] -= other[i];
    }
    return *this;
}

bool math::linal::FVector::is_zero() const {
    for (const auto& val : *this) {
        if (dcmp(val) != 0) {
            return false;
        }
    }
    return true;
}

math::linal::FVector::value_type math::linal::FVector::dot(const FVector& other) const {
    if (size() != other.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    value_type result = 0.0;
    for (size_t i = 0; i < size(); ++i) {
        result += (*this)[i] * other[i];
    }
    return result;
}

math::linal::FVector::value_type math::linal::FVector::norm() const {
    return std::sqrt(dot(*this));
}

math::linal::FVector& math::linal::FVector::normalize() {
    value_type norm = this->norm();
    if (dcmp(norm) == 0) {
        throw std::runtime_error("Cannot normalize zero vector");
    }
    return *this /= norm;
}

bool math::linal::operator==(const FVector& v1, const FVector& v2) {
    if (v1.size() != v2.size()) 
        return false;
    for (size_t i = 0; i < v1.size(); ++i) {
        if (dcmp(v1[i], v2[i]) != 0) 
            return false;
    }
    return true;
}

bool math::linal::operator!=(const FVector& v1, const FVector& v2) {
    return !(v1 == v2);
}

math::linal::FVector math::linal::operator*(FVector vector, FVector::value_type scalar) {
    return vector *= scalar;
}

math::linal::FVector math::linal::operator*(FVector::value_type scalar, FVector vector) {
    return vector *= scalar;
}

math::linal::FVector math::linal::operator/(FVector vector, FVector::value_type scalar) {
    return vector /= scalar;
}

math::linal::FVector math::linal::operator+(const FVector& v1, const FVector& v2) {
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    FVector result(v1);
    return result += v2;
}

math::linal::FVector math::linal::operator-(const FVector& v1, const FVector& v2) {
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    FVector result(v1);
    return result -= v2;
}

math::linal::FVector math::linal::operator-(const FVector& vector) {
    return -1 * vector;
}

math::linal::FVector math::linal::normalized(FVector vector) {
    FVector::value_type norm = 0.0;
    for (const auto& val : vector) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    if (dcmp(norm) == 0) {
        throw std::runtime_error("Cannot normalize zero vector");
    }
    return vector * (1.0 / norm);
}