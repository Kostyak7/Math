#pragma once 

#include <tuple>

namespace math::linal {

	template <class TMatrix>
	class SVDDecompositioner {
	public:
		std::tuple<TMatrix, TMatrix, TMatrix> decompose(const TMatrix& matrix) const;
	};

} // math::linal

#include "SVDDecompositioner.hpp"