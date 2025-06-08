file(GLOB_RECURSE PATTERNS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
if(NOT PATTERNS_SRC)
    add_library(patterns INTERFACE)
    target_include_directories(patterns INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    )
else()
    add_library(patterns ${PATTERNS_SRC})
    target_include_directories(patterns
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include
        PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../..
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    )
    generate_export_header(patterns
        BASE_NAME PATTERNS
        EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/../../patterns_export.hpp
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
