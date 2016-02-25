if (ENABLE_UMFPACK OR AMDIS_NEED_UMFPACK)
    if (WIN32)
        set(SuiteSparse_USE_LAPACK_BLAS ON)
    endif (WIN32)

    if (SuiteSparse_DIR)
       file(TO_CMAKE_PATH ${SuiteSparse_DIR} SuiteSparse_DIR)
    endif (SuiteSparse_DIR)

    # try to use a cmake-package of suitesparse
    find_package(SuiteSparse QUIET HINTS ${AMDIS_SuiteSparse_DIR})
    if (SuiteSparse_FOUND)
        message(STATUS "Found SuiteSparse CMake-library")
        include(${USE_SuiteSparse})
        target_link_libraries(amdis_base INTERFACE ${SuiteSparse_LIBRARIES})
        target_include_directories(amdis_base INTERFACE ${SuiteSparse_INCLUDE_DIR} ${SuiteSparse_METIS_INCLUDE_DIR})
        set(FOUND_SUITESPARSE_LIBS ${SuiteSparse_LIBRARIES})
    else (SuiteSparse_FOUND)
        # find umfpack manually by searching for umfpack.h header file
        find_library(UMFPACK_LIBRARY umfpack
            HINTS ${AMDIS_UMFPACK_LIB_DIR}
            DOC "Library file for UMFPACK")
        find_file(UMFPACK_H umfpack.h
            HINTS ${AMDIS_UMFPACK_INCLUDE_DIR} ENV CPATH /usr/include /usr/include/suitesparse /usr/include/ufsparse
            DOC "Headerfile umfpack.h for UMFPACK")

        if (UMFPACK_H AND UMFPACK_LIBRARY)
            get_filename_component(UMFPACK_PATH ${UMFPACK_H} PATH)
            get_filename_component(UMFPACK_LIB_PATH ${UMFPACK_LIBRARY} PATH)
            set(FOUND_SUITESPARSE_LIBS ${UMFPACK_LIBRARY})

            # find all connected libraries
            find_library(AMD_LIBRARY amd HINTS ${UMFPACK_LIB_PATH})
            find_library(BLAS_LIBRARY NAMES blas openblas HINTS ${UMFPACK_LIB_PATH} /usr/lib /usr/lib/openblas-base)
            find_library(CHOLMOD_LIBRARY cholmod HINTS ${UMFPACK_LIB_PATH})
            find_library(COLAMD_LIBRARY colamd HINTS ${UMFPACK_LIB_PATH})
            find_library(SUITESPARSECONFIG_LIBRARY suitesparseconfig HINTS ${UMFPACK_LIB_PATH})
            if (AMD_LIBRARY AND BLAS_LIBRARY)
                list(APPEND FOUND_SUITESPARSE_LIBS ${AMD_LIBRARY} ${BLAS_LIBRARY})
            endif (AMD_LIBRARY AND BLAS_LIBRARY)
            if (CHOLMOD_LIBRARY)
                list(APPEND FOUND_SUITESPARSE_LIBS ${CHOLMOD_LIBRARY})
            endif (CHOLMOD_LIBRARY)
            if (COLAMD_LIBRARY)
                list(APPEND FOUND_SUITESPARSE_LIBS ${COLAMD_LIBRARY})
            endif (COLAMD_LIBRARY)
            if (SUITESPARSECONFIG_LIBRARY)
                list(APPEND FOUND_SUITESPARSE_LIBS ${SUITESPARSECONFIG_LIBRARY})
            endif (SUITESPARSECONFIG_LIBRARY)

            target_include_directories(amdis_base INTERFACE ${UMFPACK_PATH})
            target_link_libraries(amdis_base INTERFACE ${FOUND_SUITESPARSE_LIBS})
        else()
            message(FATAL_ERROR "Could not find the UMFPACK header umfpack.h.")
        endif (UMFPACK_H AND UMFPACK_LIBRARY)
    endif (SuiteSparse_FOUND)


    # Check for clock_gettime in librt
    if (NOT WIN32)
	include(CheckLibraryExists)
	check_library_exists(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME)
	if (HAVE_CLOCK_GETTIME)
	    target_link_libraries(amdis_base INTERFACE rt)
	else ()
	    check_library_exists(c clock_gettime "" HAVE_CLOCK_GETTIME)
	endif (HAVE_CLOCK_GETTIME)
    endif (NOT WIN32)


    # collect informations about umfpack version and found libraries
    if (FOUND_SUITESPARSE_LIBS)
        find_file(_UMFPACK_H umfpack.h HINTS ${SuiteSparse_INCLUDE_DIR} ${UMFPACK_PATH})
        file(STRINGS ${_UMFPACK_H} UMFPACK_VERSION_LINE REGEX "#define UMFPACK_VERSION")
        string(REGEX MATCH "\"UMFPACK V?([^\"]+)\"" UMFPACK_VERSION_REGEX ${UMFPACK_VERSION_LINE})
        set(UMFPACK_VERSION ${CMAKE_MATCH_1})
        unset(_UMFPACK_H CACHE)

        message(STATUS "UMFPACK version: ${UMFPACK_VERSION}")
        message(STATUS "Found the following SuiteSparse libraries:")
        foreach (lib ${FOUND_SUITESPARSE_LIBS})
            message(STATUS "  ${lib}")
        endforeach ()
    endif (FOUND_SUITESPARSE_LIBS)

    target_compile_definitions(amdis_base INTERFACE
        HAVE_UMFPACK=1
        MTL_HAS_UMFPACK=1)

    install(FILES ${BASE_DIR}/cmake/find_umfpack.cmake DESTINATION share/amdis/)
endif (ENABLE_UMFPACK OR AMDIS_NEED_UMFPACK)