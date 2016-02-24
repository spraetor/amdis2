add_library(mtl4 INTERFACE)

if (IS_AMDISCONFIG)
    target_include_directories(mtl4 INTERFACE ${AMDIS_INCLUDE_DIR}/mtl4)
else ()
    target_include_directories(mtl4 INTERFACE ${BASE_DIR}/lib/mtl4)
endif (IS_AMDISCONFIG)

target_compile_features(mtl4 INTERFACE
    cxx_rvalue_references
    cxx_auto_type
    cxx_range_for
    cxx_generalized_initializers
    cxx_static_assert
    cxx_defaulted_functions)

set (CXX_ELEVEN_FEATURE_LIST "MOVE" "AUTO" "RANGEDFOR" "INITLIST" "STATICASSERT" "DEFAULTIMPL")
foreach (feature ${CXX_ELEVEN_FEATURE_LIST})
    target_compile_definitions(mtl4 INTERFACE MTL_WITH_${feature})
endforeach ()


if (ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
    target_compile_definitions(mtl4 INTERFACE MTL_WITH_OPENMP)
    target_compile_options(mtl4 INTERFACE ${OpenMP_CXX_FLAGS})
    else ()
            message(FATAL_ERROR "OpenMP not found")
    endif (OPENMP_FOUND)
endif (ENABLE_OPENMP)


if (NOT IS_AMDISCONFIG)
    # specify how to install this target:
    # -----------------------------------
    install(DIRECTORY ${BASE_DIR}/lib/mtl4/
        DESTINATION include/amdis/mtl4/
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN ".svn" EXCLUDE
        PATTERN ".svn/*" EXCLUDE
        PATTERN "mtl4/libs" EXCLUDE
        PATTERN "mtl4/extern" EXCLUDE)
endif (NOT IS_AMDISCONFIG)

target_link_libraries(amdis_base INTERFACE mtl4)
