file(GLOB_RECURSE THREADS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_library(threads ${THREADS_SRC})
target_include_directories(threads
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../..
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)
generate_export_header(threads
    BASE_NAME THREADS
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/../../threads_export.hpp
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
