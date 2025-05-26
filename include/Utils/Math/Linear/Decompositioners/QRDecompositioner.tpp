#include <Utils/Math/Linear/AnyMatrix.h>

namespace math::linal {

	template <class TMatrix>
	std::pair<TMatrix, TMatrix> QRDecompositioner::decompose(const TMatrix& matrix) const {
		const size_t n = matrix.get_width();
		TMatrix Q(n, n);
		TMatrix R(n, n);

		return { Q, R };
	}

} // math::linal
