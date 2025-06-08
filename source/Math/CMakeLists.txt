file(GLOB_RECURSE MATH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_library(math ${MATH_SRC})
target_include_directories(math
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../..
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)
generate_export_header(math
    BASE_NAME MATH
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/../../math_export.hpp
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
