template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::Matrix(std::initializer_list<std::initializer_list<value_type>> init) {
	size_t i = 0, j = 0;
	for (const auto& row : init) {
		j = 0;
		for (const auto& val : row) {
			(*this)[i][j] = val;
			++j;
		}
		++i;
	}
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator*=(value_type scalar) {
	for (auto& vec : *this) {
		vec *= scalar;
	}
	return *this;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator/=(value_type scalar) {
	return *this *= (1 / scalar);
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator+=(const Matrix& other) {
	for (size_t r = 0; r < Rows; ++r) {
		(*this)[r] += other[r];
	}
	return *this;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator-=(const Matrix& other) {
	for (size_t r = 0; r < Columns; ++r) {
		(*this)[r] += other[r];
	}
	return *this;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_square() const {
	return Rows == Columns;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_zero() const {
	for (size_t r = 0; r < Rows; ++r) {
		for (size_t c = 0; c < Columns; ++c) {
			if (dcmp((*this)[r][c]) != 0) {
				return false;
			}
		}
	}
	return true;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_identity() const {
	if (!is_square())
		return false;
	for (size_t r = 0; r < Rows; ++r) {
		for (size_t c = r; c < Columns; ++c) {
			if (r == c && dcmp((*this)[r][c], 1.) != 0)
				return false;
			if (r != c && (dcmp((*this)[r][c]) != 0 || dcmp((*this)[c][r]) != 0))
				return false;
		}
	}
	return true;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_diagonal() const {
	if (!is_square())
		return false;
	for (size_t r = 0; r < Rows; ++r) {
		for (size_t c = r + 1; c < Columns; ++c) {
			if (dcmp((*this)[r][c]) != 0 || dcmp((*this)[c][r]) != 0)
				return false;
		}
	}
	return true;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_symmetrical() const {
	if (!is_square())
		return false;
	for (size_t r = 0; r < Rows; ++r) {
		for (size_t c = r + 1; c < Columns; ++c) {
			if (dcmp((*this)[r][c], (*this)[c][r]) != 0)
				return false;
		}
	}
	return true;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_upper_triangular() const {
	return is_square() && is_trapezoidal();
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_lower_triangular() const {
	if (!is_square())
		return false;
	for (size_t r = 0; r < Rows; ++r) {
		for (size_t c = r + 1; c < Columns; ++c) {
			if (dcmp((*this)[r][c]) != 0)
				return false;
		}
	}
	return true;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_trapezoidal() const {
	for (size_t c = 0; c < Columns; ++c) {
		for (size_t r = c + 1; r < Rows; ++r) {
			if (dcmp((*this)[r][c]) != 0)
				return false;
		}
	}
	return true;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_positive_definite() const {
	if (!is_symmetrical())
		return false;

	/*auto eigenvalues = get_eigenvalues();
	return std::all_of(eigenvalues.begin(), eigenvalues.end(),
		[](const auto& pair) { return pair.second.real() > 0; });*/
	return false;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_negative_definite() const {
	if (!is_symmetrical())
		return false;

	/*auto eigenvalues = get_eigenvalues();
	return std::all_of(eigenvalues.begin(), eigenvalues.end(),
		[](const auto& pair) { return pair.second.real() < 0; });*/
	return false;
}

template <size_t Rows, size_t Columns>
size_t math::geom::Matrix<Rows, Columns>::get_width() const {
	return Columns;
}

template <size_t Rows, size_t Columns>
size_t math::geom::Matrix<Rows, Columns>::get_height() const {
	return Rows;
}

template <size_t Rows, size_t Columns>
typename math::geom::Matrix<Rows, Columns>::value_type math::geom::Matrix<Rows, Columns>::get(size_t row, size_t col) const {
	return (*this)[row][col];
}

template <size_t Rows, size_t Columns>
void math::geom::Matrix<Rows, Columns>::set(size_t row, size_t col, value_type value) {
	(*this)[row][col] = value;
}

template <size_t Rows, size_t Columns>
typename math::geom::Matrix<Rows, Columns>::value_type math::geom::Matrix<Rows, Columns>::det() const {
	if (!is_square()) {
		// throw std::runtime_error("Determinant is defined only for square matrices");
	}

	size_t n = get_width();
	if (n == 0) return 0.0;
	if (n == 1) return (*this)[0][0];
	if (n == 2) return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];

	Matrix temp(*this);
	value_type det = 1.0;

	for (size_t i = 0; i < n; ++i) {
		size_t pivot = i;
		for (size_t j = i + 1; j < n; ++j) {
			if (fabs(temp[j][i]) > fabs(temp[pivot][i])) {
				pivot = j;
			}
		}

		if (pivot != i) {
			std::swap(temp[i], temp[pivot]);
			det = -det;
		}

		if (temp[i][i] == 0)
			return 0.0;

		det *= temp[i][i];

		for (size_t j = i + 1; j < n; ++j) {
			value_type factor = temp[j][i] / temp[i][i];
			for (size_t k = i + 1; k < n; ++k) {
				temp[j][k] -= factor * temp[i][k];
			}
		}
	}

	return det;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns> math::geom::operator*(const Matrix<Rows, Columns>& matrix, typename Matrix<Rows, Columns>::value_type scalar) {
	Matrix res(matrix);
	res *= scalar;
	return res;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns> math::geom::operator*(typename Matrix<Rows, Columns>::value_type scalar, const Matrix<Rows, Columns>& matrix) {
	return matrix * scalar;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns> math::geom::operator/(const Matrix<Rows, Columns>& matrix, typename Matrix<Rows, Columns>::value_type scalar) {
	return matrix * (1. / scalar);
}

template <size_t Rows, size_t Columns>
math::geom::Vector<Rows> math::geom::operator*(const Matrix<Rows, Columns>& matrix, const Vector<Columns>& vector) {
	Vector<Rows> res{};
	for (size_t i = 0; i < Columns; ++i) {
		res[i] = matrix[i].dot(vector);
	}
	return res;
}

template <size_t Rows, size_t Columns>
math::geom::Vector<Columns> math::geom::operator*(const Vector<Rows>& vector, const Matrix<Rows, Columns>& matrix) {
	Vector<Columns> res{};
	for (size_t i = 0; i < Columns; ++i) {
		for (size_t j = 0; j < Rows; ++j) {
			res[i] += vector[j] * matrix[j][i];
		}
	}
	return res;
}

template <size_t N, size_t M, size_t K>
math::geom::Matrix<N, K> math::geom::operator*(const Matrix<N, M>& m1, const Matrix<M, K>& m2) {
	using value_type = typename Matrix<N, N>::value_type;
	Matrix<N, K> res;
	for (size_t r = 0; r < N; ++r) {
		for (size_t c = 0; c < K; ++c) {
			value_type s = 0.0;
			for (size_t i = 0; i < M; ++i) {
				s += m1[r][i] * m2[i][c];
			}
			res[r][c] = s;
		}
	}
	return res;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns> math::geom::operator+(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2) {
	Matrix<Rows, Columns> res(m1);
	res += m2;
	return res;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns> math::geom::operator-(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2) {
	Matrix<Rows, Columns> res(m1);
	res -= m2;
	return res;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns> math::geom::operator-(const Matrix<Rows, Columns>& matrix) {
	return -1 * matrix;
}

template <size_t N>
math::geom::Matrix<N, N> math::geom::inversed(const Matrix<N, N>& matrix) {
	using value_type = typename Matrix<N, N>::value_type;
	auto inverse = identity_matrix<N>();
	Matrix temp(matrix);

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
			// throw std::runtime_error("Matrix is singular and cannot be inverted");
		}

		value_type diag = temp[i][i];
		for (size_t j = 0; j < n; ++j) {
			temp[i][j] /= diag;
			inverse[i][j] /= diag;
		}

		for (size_t k = 0; k < n; ++k) {
			if (k != i && temp[k][i] != 0) {
				value_type factor = temp[k][i];
				for (size_t j = 0; j < n; ++j) {
					temp[k][j] -= factor * temp[i][j];
					inverse[k][j] -= factor * inverse[i][j];
				}
			}
		}
	}

	return inverse;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Columns, Rows> math::geom::transposed(const Matrix<Rows, Columns>& matrix) {
	Matrix<Columns, Rows> res;
	for (size_t r = 0; r < Rows; ++r) {
		for (size_t c = 0; c < Columns; ++c) {
			res[c][r] = matrix[r][c];
		}
	}
	return res;
}

template <size_t N>
math::geom::Matrix<N, N> math::geom::identity_matrix() {
	Matrix<N, N> res{};
	for (size_t i = 0; i < N; ++i) {
		res[i][i] = 1.0;
	}
	return res;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns> math::geom::elementary_matrix_unit(size_t i, size_t j, double value) {
	Matrix<N, N> res{};
	res.set(i, j, value);
	return res;
}

math::geom::Matrix2x2 math::geom::rotation_matrix_2d(double angle_rad) {
	Matrix2x2 result = identity_matrix<2>();
	double c = std::cos(angle_rad);
	double s = std::sin(angle_rad);

	result[0][0] = c;  result[0][1] = -s;
	result[1][0] = s;  result[1][1] = c;

	return result;
}

math::geom::Matrix3x3 math::geom::rotation_matrix_x(double angle_rad) {
	Matrix3x3 result = identity_matrix<3>();
	double c = std::cos(angle_rad);
	double s = std::sin(angle_rad);

	result[1][1] = c;  result[1][2] = -s;
	result[2][1] = s;  result[2][2] = c;

	return result;
}

math::geom::Matrix3x3 math::geom::rotation_matrix_y(double angle_rad) {
	Matrix3x3 result = identity_matrix<3>();
	double c = std::cos(angle_rad);
	double s = std::sin(angle_rad);

	result[1][1] = c;  result[1][2] = -s;
	result[2][1] = s;  result[2][2] = c;

	return result;
}

math::geom::Matrix3x3 math::geom::rotation_matrix_z(double angle_rad) {
	Matrix3x3 result = identity_matrix<3>();
	double c = std::cos(angle_rad);
	double s = std::sin(angle_rad);

	result[1][1] = c;  result[1][2] = -s;
	result[2][1] = s;  result[2][2] = c;

	return result;
}

math::geom::Matrix4x4 math::geom::orthographic_projection(double left, double right, double bottom, double top, double near, double far) {
	Matrix4x4 result = identity_matrix<4>();

	result[0][0] = 2.0f / (right - left);
	result[1][1] = 2.0f / (top - bottom);
	result[2][2] = -2.0f / (far - near);

	result[0][3] = -(right + left) / (right - left);
	result[1][3] = -(top + bottom) / (top - bottom);
	result[2][3] = -(far + near) / (far - near);

	return result;
}

math::geom::Matrix4x4 math::geom::perspective_projection(double fov_rad, double aspect, double near, double far) {
	Matrix4x4 result{};
	double tan_half_fov = std::tan(fov_rad / 2.0f);

	result[0][0] = 1.0f / (aspect * tan_half_fov);
	result[1][1] = 1.0f / tan_half_fov;
	result[2][2] = -(far + near) / (far - near);
	result[2][3] = -1.0f;
	result[3][2] = -(2.0f * far * near) / (far - near);

	return result;
}

math::geom::Matrix2x2 math::geom::scaling_matrix(double sx, double sy) {
	Matrix2x2 result{};
	result[0][0] = sx;
	result[1][1] = sy;
	return result;
}

math::geom::Matrix3x3 math::geom::scaling_matrix(double sx, double sy, double sz) {
	Matrix3x3 result{};
	result[0][0] = sx;
	result[1][1] = sy;
	result[2][2] = sz;
	return result;
}

math::geom::Matrix4x4 math::geom::scaling_matrix(double sx, double sy, double sz, double sw) {
	Matrix4x4 result{};
	result[0][0] = sx;
	result[1][1] = sy;
	result[2][2] = sz;
	result[3][3] = sw;
	return result;
}

math::geom::Matrix4x4 math::geom::scale_around_point(float sx, float sy, float sz, float cx, float cy, float cz) {
	Matrix4x4 translateToOrigin = translation_matrix(-cx, -cy, -cz);
	Matrix4x4 scale = scaling_matrix(sx, sy, sz, 1.0f);
	Matrix4x4 translateBack = translation_matrix(cx, cy, cz);
	return translateBack * scale * translateToOrigin;
}

math::geom::Matrix3x3 math::geom::translation_matrix(double x, double y) {
	Matrix3x3 result = identity_matrix<3>();
	result[0][2] = x;
	result[1][2] = y;
	return result;
}

math::geom::Matrix4x4 math::geom::translation_matrix(double x, double y, double z) {
	Matrix4x4 result = identity_matrix<4>();
	result[0][3] = x;
	result[1][3] = y;
	result[2][3] = z;
	return result;
}
