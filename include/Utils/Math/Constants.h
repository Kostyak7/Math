#pragma once

#include <numeric>

namespace math {

	constexpr double pi =						3.1415926535897932384626;
	constexpr double e =						2.7182818284590452353602;
	constexpr double pfi =						0.6180339887498948482045;
	constexpr double gamma_Euler_Mascheron =	0.5772156649015328606065;

	template <typename T>
	constexpr T eps = std::numeric_limits<T>::epsilon();	

} // namespace math