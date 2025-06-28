cmake_minimum_required(VERSION 3.22)
project(THREADS LANGUAGES CXX)

# Формуруем главную цель для этого CMakeLists.txt: threads
file(GLOB_RECURSE THREADS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_library(threads ${MESH2D_SOURCE})
set_target_properties(threads PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON 
    CXX_EXTENSIONS OFF)
include(GenerateExportHeader)
generate_export_header(threads)
target_include_directories(threads
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/source/)
if (BUILD_SHARED_LIBS)
    if (UNIX)
        target_compile_options(threads PRIVATE -fvisibility=hidden -fvisibility-inlines-hidden)
    elseif(WIN32)
        add_custom_command(TARGET threads POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_RUNTIME_DLLS:threads>
            $<TARGET_FILE_DIR:threads> COMMAND_EXPAND_LISTS)
    endif()
endif()

# === Тесты ===
option(THREADS_TESTING "Enable unit tests" ON)
if (THREADS_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
