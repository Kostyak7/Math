#include <Utils/Math/Linear/AnyMatrix.h>

namespace math::linal {

	template <class TMatrix>
	TMatrix CholeskyDecompositioner::decompose(const TMatrix& matrix) const {
		const size_t n = matrix.get_width();
		TMatrix L(n, n);

		return L;
	}

} // math::linal
