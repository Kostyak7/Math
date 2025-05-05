#pragma once 

namespace math::linal {

	template <class TMatrix>
	class QRDecompositioner {
	public:
		std::pair<TMatrix, TMatrix> decompose(const TMatrix& matrix) const;
	};

} // math::linal

#include "QRDecompositioner.hpp"