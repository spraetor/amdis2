project("amdis_demo")
cmake_minimum_required(VERSION 2.8)

#find_package(AMDIS REQUIRED COMPONENTS umfpack )
find_package(AMDIS REQUIRED)

if(AMDIS_FOUND)
	message(STATUS "AMDiS was found.")
	include(${AMDIS_USE_FILE})
	SET(BASIS_LIBS ${AMDIS_LIBRARIES})
endif(AMDIS_FOUND)

add_executable("ellipt" src/ellipt.cc)
target_link_libraries("ellipt" ${BASIS_LIBS})

add_executable("test1" src/test1.cc)
target_link_libraries("test1" ${BASIS_LIBS})
