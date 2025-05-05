cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(CM_TEST LANGUAGES CXX)
message(NOTICE "Your c++ compiler: ${CMAKE_CXX_COMPILER_ID}")

set(SRC_DIR ${CMAKE_SOURCE_DIR}/source)

include_directories(${SRC_DIR})

# Prepare target for compiling binary file 
# file(GLOB_RECURSE PROJECT_SOURCES ${SRC_DIR}/*.cpp)
file(GLOB_RECURSE PROJECT_SOURCES ${SRC_DIR}/Main.cpp ${SRC_DIR}/Utils/Math/Linear/*.cpp ${SRC_DIR}/Tests/*.cpp)
add_executable(${PROJECT_NAME} ${PROJECT_SOURCES} ${PROJECT_HEADERS})

# Checking of the macro __FUNCTION__ and fixing it if not supported
include(CheckCXXSourceCompiles)
check_cxx_source_compiles([[
	int main(void) {
		const char str[] = __FUNCTION__" str ";
		return 0;
	}
]] FUNCTION_MACRO_SUPPORT)
if (NOT FUNCTION_MACRO_SUPPORT)
	target_compile_definitions(${PROJECT_NAME} PUBLIC -D __FUNCTION__="")
endif()

target_include_directories(${PROJECT_NAME} PRIVATE ${SRC_DIR}/Utils)
