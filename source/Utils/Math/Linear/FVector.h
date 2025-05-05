#pragma once 

#include <vector>

namespace math::linal {

	class FVector : public std::vector<double> {
	public:
        using vector::vector;

        FVector& operator*=(value_type scalar);
        FVector& operator/=(value_type scalar);

        FVector& operator+=(const FVector& other);
        FVector& operator-=(const FVector& other);

        bool is_zero() const;

        value_type dot(const FVector& other) const;

        value_type norm() const;

        FVector& normalize();
	};

    bool operator==(const FVector& v1, const FVector& v2);
    bool operator!=(const FVector& v1, const FVector& v2);

    FVector operator*(FVector vector, FVector::value_type scalar);
    FVector operator*(FVector::value_type scalar, FVector vector);
    FVector operator/(FVector vector, FVector::value_type scalar);

    FVector operator+(const FVector& v1, const FVector& v2);
    FVector operator-(const FVector& v1, const FVector& v2);

    FVector operator-(const FVector& vector);

    FVector normalized(FVector vector);

} // namespace math::linal
