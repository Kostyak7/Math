#include <Utils/Math/Linear/AnyMatrix.h>

namespace math::linal {

	template <class TMatrix>
	std::tuple<TMatrix, TMatrix, TMatrix> SVDDecompositioner::decompose(const TMatrix& matrix) const {
		const size_t n = matrix.get_width();
		TMatrix S(n, n);
		TMatrix V(n, n);
		TMatrix D(n, n);

		return { S, V, D };
	}

} // math::linal
