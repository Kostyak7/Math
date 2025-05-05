#pragma once 

#include <iostream>

namespace tests {

	inline void print_bool(bool v) {
		std::cout << (v ? "True" : "False");
	}

} // namespace tests
