#include <Math/Linear/IMatrix.h>

#include <Math/Common.h>

#include <algorithm>
#include <random>
#include <stdexcept>

bool math::linal::IMatrix::is_equal(const IMatrix& other) const {
	if (get_width() != other.get_width() || get_height() != other.get_height())
		return false;
	for (size_t r = 0; r < get_height(); ++r) {
		for (size_t c = 0; c < get_width(); ++c) {
			if (dcmp(get(r, c), other.get(r, c)) != 0) {
				return false;
			}
		}
	}
	return true;
}

bool math::linal::IMatrix::is_empty() const {
	return get_height() == 0 || get_width() == 0;
}

bool math::linal::IMatrix::is_square() const {
	return get_width() == get_height();
}

bool math::linal::IMatrix::is_zero() const {
	for (size_t r = 0; r < get_height(); ++r) {
		for (size_t c = 0; c < get_width(); ++c) {
			if (dcmp(get(r, c)) != 0) {
				return false;
			}
		}
	}
	return true;
}

bool math::linal::IMatrix::is_identity() const {
	if (!is_square())
		return false;
	value_type tmp;
	for (size_t r = 0; r < get_height(); ++r) {
		for (size_t c = r; c < get_width(); ++c) {
			tmp = get(r, c);
			if (r == c && dcmp(get(r, c), 1.) != 0)
				return false;
			if (r != c && (dcmp(tmp) != 0 || dcmp(get(c, r)) != 0))
				return false;
		}
	}
	return true;
}

bool math::linal::IMatrix::is_diagonal() const {
	if (!is_square())
		return false;
	for (size_t r = 0; r < get_height(); ++r) {
		for (size_t c = r + 1; c < get_width(); ++c) {
			if (dcmp(get(r, c)) != 0 || dcmp(get(c, r)) != 0)
				return false;
		}
	}
	return true;
}

bool math::linal::IMatrix::is_symmetrical() const {
	if (!is_square())
		return false;
	for (size_t r = 0; r < get_height(); ++r) {
		for (size_t c = r + 1; c < get_width(); ++c) {
			if (dcmp(get(r, c), get(c, r)) != 0)
				return false;
		}
	}
	return true;
}

bool math::linal::IMatrix::is_upper_triangular() const {
	return is_square() && is_trapezoidal();
}

bool math::linal::IMatrix::is_lower_triangular() const {
	if (!is_square())
		return false;
	for (size_t r = 0; r < get_height(); ++r) {
		for (size_t c = r + 1; c < get_width(); ++c) {
			if (dcmp(get(r, c)) != 0)
				return false;
		}
	}
	return true;
}

bool math::linal::IMatrix::is_trapezoidal() const {
	for (size_t c = 0; c < get_width(); ++c) {
		for (size_t r = c + 1; r < get_height(); ++r) {
			if (dcmp(get(r, c)) != 0)
				return false;
		}
	}
	return true;
}

bool math::linal::IMatrix::is_positive_definite() const {
	if (!is_symmetrical()) 
		return false;

	auto eigenvalues = get_eigenvalues();
	return std::all_of(eigenvalues.begin(), eigenvalues.end(),
		[](const auto& pair) { return pair.second.real() > 0; });
}

bool math::linal::IMatrix::is_negative_definite() const {
	if (!is_symmetrical())
		return false;

	auto eigenvalues = get_eigenvalues();
	return std::all_of(eigenvalues.begin(), eigenvalues.end(),
		[](const auto& pair) { return pair.second.real() < 0; });
}

bool math::linal::is_positive_definite_stochastic(const IMatrix& matrix, size_t test_points) {
	if (!matrix.is_square())
		return false;

	using value_type = IMatrix::value_type;

	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::normal_distribution<value_type> dist(0, 1);

	const size_t n = matrix.get_width();
	DVector x(n);

	for (size_t test = 0; test < test_points; ++test) {
		for (size_t i = 0; i < n; ++i) {
			x[i] = dist(gen);
		}

		value_type quadratic_form = 0;
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				quadratic_form += x[i] * matrix.get(i, j) * x[j];
			}
		}

		if (quadratic_form <= 0) return false;
	}
	return true;
}

bool math::linal::is_negative_definite_stochastic(const IMatrix& matrix, size_t test_points) {
	if (!matrix.is_square())
		return false;

	using value_type = IMatrix::value_type;

	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::normal_distribution<value_type> dist(0, 1);

	const size_t n = matrix.get_width();
	DVector x(n);

	for (size_t test = 0; test < test_points; ++test) {
		for (size_t i = 0; i < n; ++i) {
			x[i] = dist(gen);
		}

		value_type quadratic_form = 0;
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				quadratic_form += x[i] * matrix.get(i, j) * x[j];
			}
		}

		if (quadratic_form >= 0) return false;
	}
	return true;
}