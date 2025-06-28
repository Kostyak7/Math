#include <Utils/Math/Linear/AnyMatrix.h>

namespace math::linal {

	template <class TMatrix>
	std::pair<TMatrix, TMatrix> LPDecompositioner::decompose(const TMatrix& matrix) const {
		const size_t n = matrix.get_width();
		TMatrix L(n, n);
		TMatrix P(n, n);

		return { L, P };
	}

} // math::linal
