cmake_minimum_required(VERSION 3.22)
project(Utils LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Включим генерацию экспортных заголовков
include(GenerateExportHeader)

# === Math ===
file(GLOB_RECURSE MATH_SRC source/Math/*.cpp)
add_library(math ${MATH_SRC})
target_include_directories(math
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/Math
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/source/Math
)
generate_export_header(math
    BASE_NAME MATH
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/math_export.hpp
)
if (BUILD_SHARED_LIBS)
    if (UNIX)
        target_compile_options(math PRIVATE -fvisibility=hidden -fvisibility-inlines-hidden)
    elseif(WIN32)
        add_custom_command(TARGET math POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_RUNTIME_DLLS:math>
            $<TARGET_FILE_DIR:math> COMMAND_EXPAND_LISTS)
    endif()
endif()


# === Patterns ===
file(GLOB_RECURSE PATTERNS_SRC source/Patterns/*.cpp)
if(NOT PATTERNS_SRC)
    add_library(patterns INTERFACE)
    target_include_directories(patterns INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Patterns
    )
else()
    add_library(patterns ${PATTERNS_SRC})
    target_include_directories(patterns
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/Patterns
        PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/source/Patterns
    )
    generate_export_header(patterns
        BASE_NAME PATTERNS
        EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/patterns_export.hpp
    )
    if (BUILD_SHARED_LIBS)
        if (UNIX)
            target_compile_options(patterns PRIVATE -fvisibility=hidden -fvisibility-inlines-hidden)
        elseif(WIN32)
            add_custom_command(TARGET patterns POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
                $<TARGET_RUNTIME_DLLS:patterns>
                $<TARGET_FILE_DIR:patterns> COMMAND_EXPAND_LISTS)
        endif()
    endif()
endif()


# === Threads ===
file(GLOB_RECURSE THREADS_SRC source/Threads/*.cpp)
add_library(threads ${THREADS_SRC})
target_include_directories(threads
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/Threads
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/source/Threads
)
generate_export_header(threads
    BASE_NAME THREADS
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/threads_export.hpp
)
if (BUILD_SHARED_LIBS)
    if (UNIX)
        target_compile_options(threads PRIVATE -fvisibility=hidden -fvisibility-inlines-hidden)
    elseif(WIN32)
        add_custom_command(TARGET threads POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_RUNTIME_DLLS:threads>
            $<TARGET_FILE_DIR:threads> COMMAND_EXPAND_LISTS)
    endif()
endif()

# === Общая interface-библиотека Utils ===
add_library(Utils INTERFACE)
target_link_libraries(Utils
  INTERFACE math patterns threads
)
target_include_directories(Utils INTERFACE
${CMAKE_CURRENT_SOURCE_DIR}/include
${CMAKE_CURRENT_BINARY_DIR}
)

# === Подключение зависимостей ===
# find_package(nlohmann_json CONFIG REQUIRED)
# target_link_libraries(Utils PUBLIC nlohmann_json::nlohmann_json) 

# === Тесты ===
option(UTILS_TESTING "Enable unit tests" ON)
if (UTILS_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()
