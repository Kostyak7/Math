#pragma once 

#include "math_export.hpp"

#include <vector>

namespace math::linal {

	class MATH_EXPORT DVector : public std::vector<double> {
	public:
        using value_type = double;
        using vector::vector;

        DVector& operator*=(value_type scalar);
        DVector& operator/=(value_type scalar);

        DVector& operator+=(const DVector& other);
        DVector& operator-=(const DVector& other);

        bool is_zero() const;

        value_type dot(const DVector& other) const;
        value_type norm() const;
        value_type length() const;

        DVector& normalize();
	};

    MATH_EXPORT bool operator==(const DVector& v1, const DVector& v2);
    MATH_EXPORT bool operator!=(const DVector& v1, const DVector& v2);

    MATH_EXPORT DVector operator*(const DVector& vector, DVector::value_type scalar);
    MATH_EXPORT DVector operator*(DVector::value_type scalar, const DVector& vector);
    MATH_EXPORT DVector operator/(const DVector& vector, DVector::value_type scalar);

    MATH_EXPORT DVector operator+(const DVector& v1, const DVector& v2);
    MATH_EXPORT DVector operator-(const DVector& v1, const DVector& v2);

    MATH_EXPORT DVector operator-(const DVector& vector);

    DVector MATH_EXPORT normalized(DVector vector);

} // namespace math::linal
