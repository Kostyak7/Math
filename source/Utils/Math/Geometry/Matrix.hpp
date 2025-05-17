//////////////////////////////////////////
//                Matrix                //
//////////////////////////////////////////

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::Matrix()
	: m_data({})
{
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::Matrix(const Matrix& matrix) {
	for (size_t i = 0; i < Rows * Columns; ++i) {
		m_data[i] = matrix.m_data[i];
	}
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::Matrix(Matrix&& matrix) noexcept 
	: m_data(std::move(matrix.m_data))
{
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::Matrix(std::initializer_list<std::initializer_list<value_type>> init) {
	// ...
	/*size_t i = 0, j = 0;
	for (const auto& row : init) {
		j = 0;
		for (const auto& val : row) {
			(*this)[i][j] = val;
			++j;
		}
		++i;
	}*/
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator=(const Matrix& matrix) {
	for (size_t i = 0; i < Rows * Columns; ++i) {
		m_data[i] = matrix.m_data[i];
	}
	return *this;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator=(Matrix&& matrix) noexcept {
	m_data = std::move(matrix.m_data);
	return *this;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator*=(value_type scalar) {
	for (auto& val : m_data) {
		val *= scalar;
	}
	return *this;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator/=(value_type scalar) {
	return *this *= (1. / scalar);
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator+=(const Matrix& other) {
	for (size_t i = 0; i < m_data.size(); ++i) {
		m_data[i] += other.m_data[i];
	}
	return *this;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>& math::geom::Matrix<Rows, Columns>::operator-=(const Matrix& other) {
	for (size_t i = 0; i < m_data.size(); ++i) {
		m_data[i] -= other.m_data[i];
	}
	return *this;
}

template <size_t Rows, size_t Columns>
void math::geom::Matrix<Rows, Columns>::swap(const Matrix& other) noexcept {
	std::swap(m_data, other.m_data);
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_empty() const {
	return Rows == 0 || Columns == 0;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_square() const {
	return Rows == Columns;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_zero() const {
	for (size_t i = 0; i < m_data.size(); ++i) {
		if (dcmp(m_data[i]) != 0)
			return false;
	}
	return true;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_identity() const {
	if (!is_square())
		return false;
	for (size_t r = 0; r < Rows; ++r) {
		for (size_t c = r; c < Columns; ++c) {
			if (r == c && dcmp(_(r, c), 1.) != 0)
				return false;
			if (r != c && (dcmp(_(r, c)) != 0 || dcmp(_(c, r)) != 0))
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
			if (dcmp(_(r, c)) != 0 || dcmp(_(c, r)) != 0)
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
			if (dcmp(_(r, c), _(c, r)) != 0)
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
			if (dcmp(_(r, c)) != 0)
				return false;
		}
	}
	return true;
}

template <size_t Rows, size_t Columns>
bool math::geom::Matrix<Rows, Columns>::is_trapezoidal() const {
	for (size_t c = 0; c < Columns; ++c) {
		for (size_t r = c + 1; r < Rows; ++r) {
			if (dcmp(_(r, c)) != 0)
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
	return _(row, col);
}

template <size_t Rows, size_t Columns>
void math::geom::Matrix<Rows, Columns>::set(size_t row, size_t col, value_type value) {
	_(row, col) = value;
}

template <size_t Rows, size_t Columns>
typename math::geom::Matrix<Rows, Columns>::ConstProxyVector math::geom::Matrix<Rows, Columns>::operator[](size_t row) const {
	return { std::next(m_data, row * Columns) };
}

template <size_t Rows, size_t Columns>
typename math::geom::Matrix<Rows, Columns>::ProxyVector math::geom::Matrix<Rows, Columns>::operator[](size_t row) {
	return { std::next(m_data, row * Columns) };
}

template <size_t Rows, size_t Columns>
typename math::geom::Matrix<Rows, Columns>::ConstProxyVector math::geom::Matrix<Rows, Columns>::at(size_t row) const {
	// ...
	return { std::next(m_data, row * Columns) };
}

template <size_t Rows, size_t Columns>
typename math::geom::Matrix<Rows, Columns>::ProxyVector math::geom::Matrix<Rows, Columns>::at(size_t row) {
	// ...
	return { std::next(m_data, row * Columns) };
}

template <size_t Rows, size_t Columns>
void math::geom::Matrix<Rows, Columns>::swap_rows(size_t r1, size_t r2) {

}

template <size_t Rows, size_t Columns>
void math::geom::Matrix<Rows, Columns>::swap_columns(size_t c1, size_t c2) {

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
void math::geom::swap(Matrix<Rows, Columns>& m1, Matrix<Rows, Columns>& m2) noexcept {
	m1.swap(m2);
}

template <size_t Rows, size_t Columns>
bool math::geom::operator==(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2) {
	if (&m1 == &m2) return true;
	for (size_t r = 0; r < Rows; ++r) {
		auto row1 = m1[r];
		for (auto it1 = row1.begin(), it2 = m2[r].begin(); it1 != row1.end(); ++it1, ++it2) {
			if (dcmp(*it1, *it2) != 0) {
				return false;
			}
		}
	}
	return true;
}
template <size_t Rows, size_t Columns>
bool math::geom::operator!=(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2) {
	return !(m1 == m2);
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
		auto row = matrix[i];
		for (size_t j = 0; j < Rows; ++j) {
			res[i] = row[j] * vector[j];
		}
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
	return -1. * matrix;
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
			temp.swap_rows(i, pivot);
			inverse.swap_rows(i, pivot);
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

//////////////////////////////////////////
//           ConstProxyVector           //
//////////////////////////////////////////

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ConstProxyVector::const_iterator math::geom::Matrix<Rows, Columns>::ConstProxyVector::begin() const {
	return m_front;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ConstProxyVector::const_iterator math::geom::Matrix<Rows, Columns>::ConstProxyVector::end() const {
	return std::next(m_front, Columns);
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ConstProxyVector::ConstProxyVector(const_iterator front)
	: m_front(front)
{
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ConstProxyVector::operator math::geom::Vector<Columns>() const {
	Vector<Columns> res;
	auto it = m_front;
	for (size_t i = 0; i < Columns; ++i, ++it) {
		res[i] = *it;
	}
	return res;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::value_type math::geom::Matrix<Rows, Columns>::ConstProxyVector::operator[](size_t col) const {
	return *std::next(m_front, col);
}

template <size_t Rows, size_t Columns>
size_t math::geom::Matrix<Rows, Columns>::ConstProxyVector::size() const {
	return Columns;
}

bool math::geom::Matrix<Rows, Columns>::ConstProxyVector::empty() const {
	return Columns == 0;
}

template <size_t Rows, size_t Columns>
bool math::geom::operator==(const Matrix<Rows, Columns>::ConstProxyVector& v1, const Matrix<Rows, Columns>::ConstProxyVector& v2) {
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

template <size_t Rows, size_t Columns>
bool math::geom::operator!=(const Matrix<Rows, Columns>::ConstProxyVector& v1, const Matrix<Rows, Columns>::ConstProxyVector& v2) {
	return !(v1 == v2);
}

//////////////////////////////////////////
//             ProxyVector              //
//////////////////////////////////////////

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ProxyVector::const_iterator math::geom::Matrix<Rows, Columns>::ProxyVector::begin() const {
	return m_front;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ProxyVector::const_iterator math::geom::Matrix<Rows, Columns>::ProxyVector::end() const {
	return std::next(m_front, Columns);
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ProxyVector::iterator math::geom::Matrix<Rows, Columns>::ProxyVector::begin() {
	return m_front;
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ProxyVector::iterator math::geom::Matrix<Rows, Columns>::ProxyVector::end() {
	return std::next(m_front, Columns);
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ProxyVector::ProxyVector(iterator front)
	: m_front(front)
{
}

template <size_t Rows, size_t Columns>
math::geom::Matrix<Rows, Columns>::ProxyVector::operator math::geom::Vector<Columns>() const {
	Vector<Columns> res;
	auto it = m_front;
	for (size_t i = 0; i < Columns; ++i, ++it) {
		res[i] = *it;
	}
	return res;
}

math::geom::Matrix<Rows, Columns>::value_type math::geom::Matrix<Rows, Columns>::ProxyVector::operator[](size_t col) const {
	return *std::next(m_front, col);
}

math::geom::Matrix<Rows, Columns>::value_type& math::geom::Matrix<Rows, Columns>::ProxyVector::operator[](size_t col) {
	return *std::next(m_front, col);
}

size_t math::geom::Matrix<Rows, Columns>::ProxyVector::size() const {
	return Columns;
}

bool math::geom::Matrix<Rows, Columns>::ProxyVector::empty() const {
	return Columns == 0;
}

bool math::geom::operator==(const Matrix<Rows, Columns>::ProxyVector& v1, const Matrix<Rows, Columns>::ProxyVector& v2) {
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

bool math::geom::operator!=(const Matrix<Rows, Columns>::ProxyVector& v1, const Matrix<Rows, Columns>::ProxyVector& v2) {
	return !(v1 == v2);
}