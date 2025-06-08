cmake_minimum_required(VERSION 3.22)
project(Utils LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Включим генерацию экспортных заголовков
include(GenerateExportHeader)

# Опции для сборки подмодулей
option(BUILD_UTILS_MATH "Build Math module" ON)
option(BUILD_UTILS_PATTERNS "Build Patterns module" ON)
option(BUILD_UTILS_THREADS "Build Threads module" ON)

# Добавляем подмодули
if(BUILD_UTILS_MATH)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/source/Math)
endif()

if(BUILD_UTILS_PATTERNS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/source/Patterns)
endif()

if(BUILD_UTILS_THREADS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/source/Threads)
endif()

# === Общая interface-библиотека Utils ===
add_library(Utils INTERFACE)
target_include_directories(Utils INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Подключаем только те модули, которые были собраны
if(BUILD_UTILS_MATH)
    target_link_libraries(Utils INTERFACE math)
endif()
if(BUILD_UTILS_PATTERNS)
    target_link_libraries(Utils INTERFACE patterns)
endif()
if(BUILD_UTILS_THREADS)
    target_link_libraries(Utils INTERFACE threads)
endif()

# === Тесты ===
option(UTILS_TESTING "Enable unit tests" ON)
if (UTILS_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
