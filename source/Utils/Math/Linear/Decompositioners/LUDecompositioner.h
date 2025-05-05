#pragma once 

namespace math::linal {

	template <class TMatrix>
	class LUDecompositioner {
	public:
		std::pair<TMatrix, TMatrix> decompose(const TMatrix& matrix) const;
	};

} // math::linal

#include "LUDecompositioner.hpp"