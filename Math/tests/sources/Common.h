#pragma once 

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <fstream>
#include <thread>
#include <string>
#include <chrono>
#include <cstdlib>
#include <random>
#include <limits>
#include <iostream>

namespace tests {

	inline void print_bool(bool v) {
		std::cout << (v ? "True" : "False");
	}

} // namespace tests
