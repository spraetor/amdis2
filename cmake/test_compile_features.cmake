
macro(test_compile_features RESULT_VAR TARGET SCOPE) # feature1 feature2 feature3 ...)
    set(_FEATURES "")
    foreach(FEATURE ${ARGN})
        set(_FEATURES "${_FEATURES} ${FEATURE}")
    endforeach()
    set(_DIR ${CMAKE_BINARY_DIR}/test_cxx11)

    # test only of variable RESULT_VAR not yet set manually or by a
    # previous call to this macro
    if (NOT DEFINED ${RESULT_VAR})
        message(STATUS "Test for cxx11 features for target ${TARGET}")
        execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${_DIR})
        file(WRITE ${_DIR}/CMakeLists.txt "
            cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
            project(TestCxx LANGUAGES CXX)
            add_library(test_cxx ${BASE_DIR}/cmake3/dummy.cc)
            target_compile_features(test_cxx PRIVATE ${_FEATURES})")

        execute_process(COMMAND ${CMAKE_COMMAND} ${_DIR}
            WORKING_DIRECTORY ${_DIR}
            RESULT_VARIABLE _ERR_CODE OUTPUT_QUIET)

        execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${_DIR})
    elseif (${RESULT_VAR})
        set(_ERR_CODE "0")
    else ()
        set(_ERR_CODE "1")
    endif ()

    if ("${_ERR_CODE}" STREQUAL "0")
        set(${RESULT_VAR} true CACHE BOOL "Enable C++11 compiler features")
        target_compile_features(${TARGET} ${SCOPE} ${ARGN})
    else ()
        set(${RESULT_VAR} false CACHE BOOL "Enable C++11 compiler features")
    endif ()
    mark_as_advanced(${RESULT_VAR})

    unset(_ERR_CODE)
    unset(_DIR)
    unset(_FEATURES)
endmacro(test_compile_features)