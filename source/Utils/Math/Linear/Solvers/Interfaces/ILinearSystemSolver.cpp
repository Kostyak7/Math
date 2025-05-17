#include "ILinearSystemSolver.h"

math::linal::ILinearSystemSolver::ILinearSystemSolver(const Params& params)
	: m_params(params)
{
}

math::linal::ILinearSystemSolver::Params math::linal::ILinearSystemSolver::get_params() const {
	return m_params;
}

bool math::linal::ILinearSystemSolver::slae_check(const IMatrix& matrix, const DVector& rhs) const {
	if (matrix.is_empty()) {
		if (m_params.throw_exceptions)
			throw std::invalid_argument("Matrix of system is empty");
		return false;
	}
	if (!matrix.is_square()) {
		if (m_params.throw_exceptions)
			throw std::invalid_argument("Matrix of system is not square");
		return false;
	}
	if (matrix.get_height() != rhs.size()) {
		if (m_params.throw_exceptions)
			throw std::invalid_argument("Matrix and rhs dimensions do not match");
		return false;
	}

	return true;
}