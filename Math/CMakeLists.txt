cmake_minimum_required(VERSION 3.22)
project(MATH LANGUAGES CXX)

# Формуруем главную цель для этого CMakeLists.txt: math
file(GLOB_RECURSE MATH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_library(math ${MESH2D_SOURCE})
set_target_properties(math PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON 
    CXX_EXTENSIONS OFF)
include(GenerateExportHeader)
generate_export_header(math)
target_include_directories(math
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/source/)
if (BUILD_SHARED_LIBS)
    if (UNIX)
        target_compile_options(math PRIVATE -fvisibility=hidden -fvisibility-inlines-hidden)
    elseif(WIN32)
        add_custom_command(TARGET math POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_RUNTIME_DLLS:math>
            $<TARGET_FILE_DIR:math> COMMAND_EXPAND_LISTS)
    endif()
endif()

# === Тесты ===
option(MATH_TESTING "Enable unit tests" ON)
if (MATH_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
