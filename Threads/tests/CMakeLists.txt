cmake_minimum_required(VERSION 3.22)

project(utils_tests LANGUAGES CXX)
find_package(GTest REQUIRED CONFIG)

# Подключаем наш CMake модуль для удобства дальнейшей линковки тестов:
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules)
include("Test")

# Перечисление всех категорий тестов:
test(utils_math_test_linear ${CMAKE_CURRENT_SOURCE_DIR}/sources/Math/Linear/*.cpp)
test(utils_math_test_geometry ${CMAKE_CURRENT_SOURCE_DIR}/sources/Math/Geometry/*.cpp)

test(utils_patterns_test_patterns ${CMAKE_CURRENT_SOURCE_DIR}/sources/Patterns/*.cpp)

test(utils_threads_test_threads ${CMAKE_CURRENT_SOURCE_DIR}/sources/Threads/*.cpp)

# Запуск всех тестов разом:
test(utils_ALL_tests ${CMAKE_CURRENT_SOURCE_DIR}/sources/*.cpp)