project("example")
cmake_minimum_required(VERSION 2.8)

find_package(AMDIS REQUIRED)
if(AMDIS_FOUND)
	message(STATUS "AMDiS was found.")
	include(${AMDIS_USE_FILE})
	SET(BASIS_LIBS ${AMDIS_LIBRARIES})
endif(AMDIS_FOUND)

add_executable("example" example.cc)
target_link_libraries("example" ${BASIS_LIBS})

