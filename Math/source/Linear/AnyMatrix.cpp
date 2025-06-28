#include <Math/Linear/AnyMatrix.h>

#include <stdexcept>

math::linal::AnyMatrix math::linal::operator*(const AnyMatrix& matrix, IMatrix::value_type scalar) {
	return std::visit([&scalar](auto matrix) -> AnyMatrix {
		matrix *= scalar;
		return matrix;
		}, matrix);
}

math::linal::AnyMatrix math::linal::operator*(IMatrix::value_type scalar, const AnyMatrix& matrix) {
	return matrix * scalar;
}

math::linal::AnyMatrix math::linal::operator/(const AnyMatrix& matrix, IMatrix::value_type scalar) {
	return matrix * (1 / scalar);
}

math::linal::DVector math::linal::operator*(const AnyMatrix& matrix, const DVector& vector) {
	return std::visit([&vector](const auto& matrix) -> DVector {
		if (matrix.get_width() != vector.size())
			throw std::invalid_argument("The matrix and the vector must be the same size");

		return matrix * vector;
		}, matrix);
}

math::linal::DVector math::linal::operator*(const DVector& vector, const AnyMatrix& matrix) {
	return std::visit([&vector](const auto& matrix) -> DVector {
		if (matrix.get_height() != vector.size())
			throw std::invalid_argument("The matrix and the vector must be the same size");

		return vector * matrix;
		}, matrix);
}

math::linal::AnyMatrix math::linal::operator*(const AnyMatrix& m1, const AnyMatrix& m2) {
	return std::visit([](const auto& m1, const auto& m2) -> AnyMatrix {
		if constexpr (!std::is_same_v<decltype(m1), decltype(m2)>)
			throw std::invalid_argument("Cannot multiply different matrix types");
		if (m1.get_width() != m2.get_height()) 
			throw std::invalid_argument("The matrices must be the same size");
		
		return m1 * m2;
		}, m1, m2);
}

math::linal::AnyMatrix math::linal::operator+(const AnyMatrix& m1, const AnyMatrix& m2) {
	return std::visit([](const auto& m1, const auto& m2) -> AnyMatrix {
		if constexpr (!std::is_same_v<decltype(m1), decltype(m2)>)
			throw std::invalid_argument("Cannot add different matrix types");
		if (m1.get_width() != m2.get_width() || m1.get_height() != m2.get_height())
			throw std::invalid_argument("The matrices must be the same size");
		
		return m1 + m2;
		}, m1, m2);
}

math::linal::AnyMatrix math::linal::operator-(const AnyMatrix& m1, const AnyMatrix& m2) {
	return std::visit([](const auto& m1, const auto& m2) -> AnyMatrix {
		if constexpr (!std::is_same_v<decltype(m1), decltype(m2)>)
			throw std::invalid_argument("Cannot subtract different matrix types");
		if (m1.get_width() != m2.get_width() || m1.get_height() != m2.get_height())
			throw std::invalid_argument("The matrices must be the same size");
		
		return m1 - m2;
		}, m1, m2);
}