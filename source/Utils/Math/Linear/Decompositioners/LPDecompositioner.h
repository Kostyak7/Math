#pragma once 

namespace math::linal {

	template <class TMatrix>
	class LPDecompositioner {
	public:
		std::pair<TMatrix, TMatrix> decompose(const TMatrix& matrix) const;
	};

} // math::linal

#include "LPDecompositioner.hpp"