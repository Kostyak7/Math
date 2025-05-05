#include <Utils/Math/Linear/AnyMatrix.h>

namespace math::linal {

	template <class TMatrix>
	std::pair<TMatrix, TMatrix> LUDecompositioner::decompose(const TMatrix& matrix) const {
		const size_t n = matrix.get_width();
		TMatrix L(n, n);
		TMatrix U(n, n);

		return { L, U };
	}

} // math::linal
