function(test TARGET SOURCES)
    file(GLOB_RECURSE TEST_SOURCES ${SOURCES})
    add_executable(${TARGET} ${TEST_SOURCES})
    set_target_properties(${TARGET} PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON 
        CXX_EXTENSIONS OFF)
    target_link_libraries(${TARGET} PUBLIC utils GTest::gtest_main)
    target_compile_definitions(${TARGET} PRIVATE ASSETS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/")

    #include(GoogleTest) # TODO: зачем это?
    #gtest_discover_tests(${TARGET}) # TODO: зачем это?

    if(BUILD_SHARED_LIBS)
		add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy 
								$<TARGET_RUNTIME_DLLS:${TARGET}> 
								$<TARGET_FILE_DIR:${TARGET}> COMMAND_EXPAND_LISTS)
	endif()
endfunction()
