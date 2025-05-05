#pragma once 

namespace math::linal {

	template <class TMatrix>
	class CholeskyDecompositioner {
	public:
		TMatrix decompose(const TMatrix& matrix) const;
	};

} // math::linal

#include "CholeskyDecompositioner.hpp"