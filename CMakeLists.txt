project(AMDIS)

cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0017 OLD)
if(POLICY CMP0042)
	cmake_policy(SET CMP0042 NEW)
endif()
	
SET(LIB_DIR ${AMDIS_SOURCE_DIR}/lib)
SET(SOURCE_DIR ${AMDIS_SOURCE_DIR}/src)

#TODO: use the cmake build type
SET(MTL_INCLUDE_DIR ${LIB_DIR}/mtl4/ CACHE PATH "mtl4 directory")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DDEBUG=0") #-Wno-unused-but-set-variable
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG=1 -O0") #-Wno-unused-but-set-variable

#define the build type, empty can be everything and nothing
if(CMAKE_BUILD_TYPE STREQUAL "")
	set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_CXX_COMPILER MATCHES ".*icpc")
	Message("Compiler: Intel")
	list(APPEND CMAKE_CXX_FLAGS "-diag-disable 654 -diag-disable 858") 
endif()

list(APPEND COMPILEFLAGS "-std=c++11 -Wall -pedantic")

SET(ENABLE_PARALLEL_DOMAIN "OFF" CACHE STRING "use parallel domain decomposition. please set to one of: PMTL, PETSC, OFF" )
option(ENABLE_ZOLTAN false)
option(ENABLE_SEQ_PETSC "Use sequential PETSc solvers" false)
option(ENABLE_UMFPACK "Use of UMFPACK solver" false)
option(ENABLE_HYPRE "Use HYPRE AMG solver" false)
option(ENABLE_COMPRESSION "Use output compression" false)

find_package(Boost)
if(Boost_FOUND)
  list(APPEND AMDIS_INCLUDE_DIRS ${Boost_INCLUDE_DIR})
  message("boost lib-dirs: ${Boost_LIBRARY_DIRS}")
  message("use the boost dir: ${Boost_INCLUDE_DIR}")
  if(WIN32)
	message("the windows find_boost does not set the boost library paths, please set it")
	SET(Boost_LIBRARY_DIRS CACHE PATH "The directory containing the boost libraries")
  endif(WIN32)
  list(APPEND AMDIS_LINK_DIRS ${Boost_LIBRARY_DIRS})
endif(Boost_FOUND)


SET(AMDIS_SRC ${SOURCE_DIR}/AdaptBase.cc
	      ${SOURCE_DIR}/AdaptInfo.cc
	      ${SOURCE_DIR}/AdaptInstationary.cc
	      ${SOURCE_DIR}/AdaptStationary.cc
	      ${SOURCE_DIR}/AMDiS.cc
	      ${SOURCE_DIR}/Assembler.cc
	      ${SOURCE_DIR}/BasisFunction.cc
	      ${SOURCE_DIR}/Boundary.cc
	      ${SOURCE_DIR}/BoundaryManager.cc
	      ${SOURCE_DIR}/BoundaryObject.cc
	      ${SOURCE_DIR}/CoarseningManager.cc
	      ${SOURCE_DIR}/CoarseningManager1d.cc
	      ${SOURCE_DIR}/CoarseningManager2d.cc
	      ${SOURCE_DIR}/CoarseningManager3d.cc
	      ${SOURCE_DIR}/ComponentTraverseInfo.cc
	      ${SOURCE_DIR}/CouplingIterationInterface.cc
	      ${SOURCE_DIR}/CreatorMap.cc
	      ${SOURCE_DIR}/DOFAdmin.cc
#	      ${SOURCE_DIR}/DOFIndexed.cc
	      ${SOURCE_DIR}/DOFMatrix.cc
	      ${SOURCE_DIR}/DOFVector.cc
#	      ${SOURCE_DIR}/Debug.cc
	      ${SOURCE_DIR}/DirichletBC.cc
	      ${SOURCE_DIR}/DualTraverse.cc
	      ${SOURCE_DIR}/ElInfo.cc
	      ${SOURCE_DIR}/ElInfo1d.cc
	      ${SOURCE_DIR}/ElInfo2d.cc
	      ${SOURCE_DIR}/ElInfo3d.cc
	      ${SOURCE_DIR}/ElInfoStack.cc
	      ${SOURCE_DIR}/Element.cc
	      ${SOURCE_DIR}/ElementData.cc
	      ${SOURCE_DIR}/ElementDofIterator.cc
	      ${SOURCE_DIR}/ElementRegion_ED.cc	
	      ${SOURCE_DIR}/FiniteElemSpace.cc
	      ${SOURCE_DIR}/FirstOrderAssembler.cc
	      ${SOURCE_DIR}/FirstOrderTerm.cc
	      ${SOURCE_DIR}/FixVec.cc
	      ${SOURCE_DIR}/Global.cc
	      ${SOURCE_DIR}/Initfile.cc
	      ${SOURCE_DIR}/Lagrange.cc
	      ${SOURCE_DIR}/LeafData.cc
	      ${SOURCE_DIR}/Line.cc
	      ${SOURCE_DIR}/Log.cc
	      ${SOURCE_DIR}/MacroElement.cc
	      ${SOURCE_DIR}/Marker.cc
	      ${SOURCE_DIR}/Mesh.cc
	      ${SOURCE_DIR}/MeshStructure.cc
	      ${SOURCE_DIR}/Operator.cc
	      ${SOURCE_DIR}/OperatorTerm.cc
	      ${SOURCE_DIR}/Parametric.cc
	      ${SOURCE_DIR}/PeriodicBC.cc
	      ${SOURCE_DIR}/ProblemInstat.cc
	      ${SOURCE_DIR}/ProblemStat.cc
	      ${SOURCE_DIR}/ProblemStatDbg.cc
	      ${SOURCE_DIR}/Projection.cc
	      ${SOURCE_DIR}/QPInfo.cc
	      ${SOURCE_DIR}/QPsiPhi.cc
	      ${SOURCE_DIR}/Quadrature.cc
	      ${SOURCE_DIR}/RCNeighbourList.cc
	      ${SOURCE_DIR}/RefinementManager.cc
	      ${SOURCE_DIR}/RefinementManager1d.cc
	      ${SOURCE_DIR}/RefinementManager2d.cc
	      ${SOURCE_DIR}/RefinementManager3d.cc
	      ${SOURCE_DIR}/RobinBC.cc
	      ${SOURCE_DIR}/ScalableQuadrature.cc
	      ${SOURCE_DIR}/SecondOrderAssembler.cc
	      ${SOURCE_DIR}/SecondOrderTerm.cc
	      ${SOURCE_DIR}/StandardProblemIteration.cc
	      ${SOURCE_DIR}/SubAssembler.cc
	      ${SOURCE_DIR}/SubElInfo.cc
	      ${SOURCE_DIR}/SubQuadrature.cc
	      ${SOURCE_DIR}/SurfaceQuadrature.cc
	      ${SOURCE_DIR}/SurfaceRegion_ED.cc
	      ${SOURCE_DIR}/SystemVector.cc
	      ${SOURCE_DIR}/Tetrahedron.cc
	      ${SOURCE_DIR}/Timer.cc
	      ${SOURCE_DIR}/Traverse.cc
	      ${SOURCE_DIR}/Triangle.cc
	      ${SOURCE_DIR}/VertexVector.cc
	      ${SOURCE_DIR}/ZeroOrderAssembler.cc
#	      ${SOURCE_DIR}/ZeroOrderTerm.cc
	      ${SOURCE_DIR}/est/Estimator.cc
#	      ${SOURCE_DIR}/est/RecoveryEstimator.cc
	      ${SOURCE_DIR}/est/ResidualEstimator.cc
	      ${SOURCE_DIR}/est/SimpleResidualEstimator.cc
	      ${SOURCE_DIR}/io/ArhReader.cc
	      ${SOURCE_DIR}/io/detail/ArhReader.cc
          ${SOURCE_DIR}/io/Arh2Reader.cc
          ${SOURCE_DIR}/io/Arh2Writer.cc
          ${SOURCE_DIR}/io/Arh3Reader.cc
          ${SOURCE_DIR}/io/Arh3Writer.cc
          ${SOURCE_DIR}/io/detail/Arh2Reader.cc
	      ${SOURCE_DIR}/io/detail/ArhWriter.cc
          ${SOURCE_DIR}/io/detail/Arh2Writer.cc
          ${SOURCE_DIR}/io/detail/Arh3Reader.cc
          ${SOURCE_DIR}/io/detail/Arh3Writer.cc
	      ${SOURCE_DIR}/io/DofWriter.cc
	      ${SOURCE_DIR}/io/ElementFileWriter.cc
	      ${SOURCE_DIR}/io/FileWriterInterface.cc
	      ${SOURCE_DIR}/io/FileWriter.cc
	      ${SOURCE_DIR}/io/GNUPlotWriter.cc
	      ${SOURCE_DIR}/io/MacroInfo.cc
	      ${SOURCE_DIR}/io/MacroReader.cc
	      ${SOURCE_DIR}/io/MacroWriter.cc
	      ${SOURCE_DIR}/io/PngReader.cc
	      ${SOURCE_DIR}/io/PngWriter.cc
	      ${SOURCE_DIR}/io/PovrayWriter.cc
	      ${SOURCE_DIR}/io/Spreadsheet.cc
	      ${SOURCE_DIR}/io/ValueReader.cc
	      ${SOURCE_DIR}/io/ValueWriter.cc
#	      ${SOURCE_DIR}/io/VtkWriter.cc
#	      ${SOURCE_DIR}/io/VtkVectorWriter.cc
#         ${SOURCE_DIR}/io/detail/VtkWriter.cc
#	      ${SOURCE_DIR}/nonlin/ProblemNonLin.cc
	      ${SOURCE_DIR}/solver/SolverMatrix.cc
	      ${SOURCE_DIR}/solver/Mapper.cc
#	      ${SOURCE_DIR}/time/RosenbrockAdaptInstationary.cc
#	      ${SOURCE_DIR}/time/RosenbrockMethod.cc
#	      ${SOURCE_DIR}/time/RosenbrockStationary.cc
	      )

if(ENABLE_PARALLEL_DOMAIN)
	
	### Search for MPI compiler and get compile flags and include paths. ###

	find_package(MPI REQUIRED)
	if(MPI_FOUND)
	  list(APPEND COMPILEFLAGS "${MPI_COMPILE_FLAGS}")
	  include_directories(${MPI_INCLUDE_PATH})
	endif(MPI_FOUND)
 


	### Search for ParMetis ###

	find_file(PARMETIS_HEADER_FILE "parmetis.h" PATH ENV CPATH)
	if(PARMETIS_HEADER_FILE)
		get_filename_component(PARMETIS_INCLUDE_PATH "${PARMETIS_HEADER_FILE}" PATH CACHE)
		include_directories(${PARMETIS_INCLUDE_PATH})
	else()
		message(FATAL_ERROR "Could not find ParMetis header file!")
	endif(PARMETIS_HEADER_FILE)

	find_file(PARMETIS_LIB_FILE "libparmetis.a" PATH ENV LIBRARY_PATH)
	if(PARMETIS_LIB_FILE)
		get_filename_component(PARMETIS_LINK_PATH "${PARMETIS_LIB_FILE}" PATH CACHE)
		link_directories(${PARMETIS_LINK_PATH})
	else()
		message(FATAL_ERROR "Could not find ParMetis library file!")
	endif(PARMETIS_LIB_FILE)


	
	### Search for Zoltan ###

	if(ENABLE_ZOLTAN)
	  find_file(ZOLTAN_HEADER_FILE "zoltan_cpp.h" PATH ENV CPATH)
	  if(ZOLTAN_HEADER_FILE)
		get_filename_component(ZOLTAN_HEADER_DIR "${ZOLTAN_HEADER_FILE}" PATH CACHE)
		include_directories(${ZOLTAN_HEADER_DIR})
	  else()
		message(FATAL_ERROR "Could not find Zoltan include file!")
	  endif(ZOLTAN_HEADER_FILE)

	  list(APPEND COMPILEFLAGS "-DHAVE_ZOLTAN")
	endif(ENABLE_ZOLTAN)



	list(APPEND COMPILEFLAGS "-DHAVE_PARALLEL_DOMAIN_AMDIS=1")
	SET(PARALLEL_DOMAIN_AMDIS_SRC
		${SOURCE_DIR}/parallel/BddcMlSolver.cc
               	${SOURCE_DIR}/parallel/DofComm.cc
		${SOURCE_DIR}/parallel/CheckerPartitioner.cc
		${SOURCE_DIR}/parallel/ElementObjectDatabase.cc
		${SOURCE_DIR}/parallel/InteriorBoundary.cc
		${SOURCE_DIR}/parallel/MeshDistributor.cc 
		${SOURCE_DIR}/parallel/MeshLevelData.cc 
		${SOURCE_DIR}/parallel/MeshManipulation.cc
		${SOURCE_DIR}/parallel/MeshPartitioner.cc
		${SOURCE_DIR}/parallel/MpiHelper.cc
		${SOURCE_DIR}/parallel/ParallelDebug.cc
		${SOURCE_DIR}/parallel/ParallelDofMapping.cc
		${SOURCE_DIR}/parallel/ParallelProblemStat.cc
		${SOURCE_DIR}/parallel/ParallelSolver.cc
		${SOURCE_DIR}/parallel/PeriodicMap.cc
		${SOURCE_DIR}/parallel/ParMetisPartitioner.cc
		${SOURCE_DIR}/parallel/StdMpi.cc
		${SOURCE_DIR}/parallel/ZoltanPartitioner.cc)
		
	if(ENABLE_PARALLEL_DOMAIN STREQUAL "PETSC")
		set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};/usr/share/cmake-2.8/Modules/;${CMAKE_SOURCE_DIR}/")
		set(PETSC_EXECUTABLE_RUNS ON)
		find_package(PETSc REQUIRED)
		include_directories(${PETSC_DIR}/include ${PETSC_DIR}/${PETSC_ARCH}/include)
		list(APPEND COMPILEFLAGS "-DHAVE_PARALLEL_PETSC" "-DPETSC_VERSION=${PETSC_VERSION}")
		list(APPEND AMDIS_INCLUDE_DIRS ${PETSC_DIR}/include ${PETSC_DIR}/${PETSC_ARCH}/include)
		list(APPEND PARALLEL_DOMAIN_AMDIS_SRC 
			${SOURCE_DIR}/parallel/MatrixNnzStructure.cc
			${SOURCE_DIR}/parallel/ParallelCoarseSpaceSolver.cc
			${SOURCE_DIR}/parallel/PetscHelper.cc
			${SOURCE_DIR}/parallel/PetscSolver.cc
			${SOURCE_DIR}/parallel/PetscSolverFeti.cc
			${SOURCE_DIR}/parallel/PetscSolverFetiDebug.cc
			${SOURCE_DIR}/parallel/PetscSolverFetiMonitor.cc
			${SOURCE_DIR}/parallel/PetscSolverFetiOperators.cc
			${SOURCE_DIR}/parallel/PetscSolverFetiTimings.cc
			${SOURCE_DIR}/parallel/PetscSolverGlobalMatrix.cc
			${SOURCE_DIR}/parallel/PetscSolverGlobalBlockMatrix.cc
			${SOURCE_DIR}/parallel/PetscSolverNavierStokes.cc
			${SOURCE_DIR}/parallel/PetscSolverNSCH.cc
			${SOURCE_DIR}/parallel/PetscSolverCahnHilliard2.cc
			${SOURCE_DIR}/parallel/PetscSolverCahnHilliard.cc                        
			${SOURCE_DIR}/parallel/PetscSolverSchur.cc
			${SOURCE_DIR}/solver/PetscTypes.cc)
	elseif(ENABLE_PARALLEL_DOMAIN STREQUAL "PMTL")
		set(MTL_INCLUDE_DIR "")
		find_package(MTL REQUIRED)
		list(APPEND COMPILEFLAGS "-DMTL_HAS_MPI" "-DHAVE_PARALLEL_MTL4")
# 		list(APPEND PARALLEL_DOMAIN_AMDIS_SRC 
# 			${SOURCE_DIR}/parallel/PMTL_ProblemStat.cc
# 			)
		list(APPEND AMDIS_LIBS ${MTL_LIBRARIES})
	else()
		message(SEND_ERROR "Set the domain-decomposition version correctly.")
	endif()

	INSTALL(FILES ResolveCompilerPaths.cmake FindPackageMultipass.cmake FindPETSc.cmake
		DESTINATION share/amdis/)
	if(USE_PETSC_DEV)
	  list(APPEND COMPILEFLAGS "-DHAVE_PETSC_DEV")
	endif( USE_PETSC_DEV)
endif(ENABLE_PARALLEL_DOMAIN)


if(ENABLE_SEQ_PETSC)
  find_package(MPI REQUIRED)
  if(MPI_FOUND)
    list(APPEND COMPILEFLAGS "${MPI_COMPILE_FLAGS}")
    include_directories(${MPI_INCLUDE_PATH})
  endif(MPI_FOUND)
	
  set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};/usr/share/cmake-2.8/Modules/;${CMAKE_SOURCE_DIR}/")
  set(PETSC_EXECUTABLE_RUNS ON)
  find_package(PETSc REQUIRED)
  include_directories(${PETSC_DIR}/include ${PETSC_DIR}/${PETSC_ARCH}/include)
  list(APPEND AMDIS_INCLUDE_DIRS ${PETSC_DIR}/include ${PETSC_DIR}/${PETSC_ARCH}/include)
  list(APPEND COMPILEFLAGS "-DHAVE_SEQ_PETSC" "-DPETSC_VERSION=${PETSC_VERSION}")
  
  list(APPEND AMDIS_SRC   
	${SOURCE_DIR}/solver/PetscSolver.cc
	${SOURCE_DIR}/solver/PetscTypes.cc)
endif(ENABLE_SEQ_PETSC)

if(ENABLE_UMFPACK)
	find_file(UMFPACK_H umfpack.h
			    HINTS ENV CPATH /usr/include /usr/include/suitesparse /usr/include/ufsparse
			    DOC "headerfile umfpack.h for UMFPACK")
	if(UMFPACK_H)
		get_filename_component(UMFPACK_PATH ${UMFPACK_H} PATH)
		include_directories(${UMFPACK_PATH} ${UFCONFIG_PATH} ${AMD_PATH})
		list(APPEND COMPILEFLAGS "-DHAVE_UMFPACK=1" "-DMTL_HAS_UMFPACK")
	else()
		message(FATAL_ERROR "Could not find the UMFPACK header umfpack.h.")
	endif()
	SET(RPM_DEPEND_STR "blas")
endif(ENABLE_UMFPACK)


if(ENABLE_HYPRE)
  include(HYPREConfig.cmake)
  message("have hypre: ${HAVE_HYPRE}")
  if (HAVE_HYPRE)
    if (NOT MPI_FOUND)
      find_package(MPI REQUIRED)
      if (MPI_FOUND)
	list(APPEND COMPILEFLAGS "${MPI_COMPILE_FLAGS}")
	include_directories(${MPI_INCLUDE_PATH})
      endif()
    endif()
    include_directories(${HYPRE_INCLUDE_DIRECTORIES})
    list(APPEND COMPILEFLAGS "-DMTL_HAS_HYPRE")
    list(APPEND MYLIBS ${HYPRE_LIBRARIES} )
  endif()
endif(ENABLE_HYPRE)
  
  
if(ENABLE_PNG)
	find_file(PNG_H png.h
			    HINTS ENV CPATH /usr/include
			    DOC "headerfile png.h for PNG-READER")
	if(PNG_H)
		get_filename_component(PNG_PATH ${PNG_H} PATH)
		include_directories(${PNG_PATH})
		list(APPEND COMPILEFLAGS "-DHAVE_PNG=1")
	else()
		message(FATAL_ERROR "Could not find the PNG header png.h.")
	endif()
endif(ENABLE_PNG)

#if(ENABLE_COMPRESSION)
#  list(APPEND COMPILEFLAGS "-DHAVE_COMPRESSION=1")
#endif(ENABLE_COMPRESSION)

#if(NOT ENABLE_OUTPUT)
#       message(WARNING "AMDiS cout output disabled!")
#       list(APPEND COMPILEFLAGS "-DSUPPRESS_OUTPUT")
#endif(NOT ENABLE_OUTPUT)



if(ENABLE_OPENMP)
  if(CMAKE_CXX_COMPILER MATCHES ".*icpc")
	list(APPEND COMPILEFLAGS "-openmp")
  else()
        list(APPEND COMPILEFLAGS "-fopenmp")
  endif()
endif(ENABLE_OPENMP)


#muparser includes
set(MUPARSER_SOURCE_DIR ${AMDIS_SOURCE_DIR}/lib/muparser_v2_2_4)
file(GLOB MUPARSER_SRC ${MUPARSER_SOURCE_DIR}/src/*.cpp)
list(APPEND AMDIS_INCLUDE_DIRS ${MUPARSER_SOURCE_DIR}/include)

#mtl4 includes
list(APPEND AMDIS_INCLUDE_DIRS ${MTL_INCLUDE_DIR})
list(APPEND AMDIS_INCLUDE_DIRS ${SOURCE_DIR})

include_directories(${AMDIS_INCLUDE_DIRS})

if(WIN32)
  add_library(amdis STATIC ${AMDIS_SRC} ${PARALLEL_DOMAIN_AMDIS_SRC})
  add_library(muparser STATIC ${MUPARSER_SRC})
else()
  add_library(amdis SHARED ${AMDIS_SRC} ${PARALLEL_DOMAIN_AMDIS_SRC})
  add_library(muparser SHARED ${MUPARSER_SRC})
endif()

list(APPEND AMDIS_LIBS amdis ${Boost_LIBRARIES})

if(WIN32)
	list(APPEND COMPILEFLAGS "-D_SCL_SECURE_NO_WARNINGS" "-D_CRT_SECURE_NO_WARNINGS")
endif(WIN32)


if(CMAKE_SYSTEM_NAME MATCHES Darwin)
  if(CMAKE_COMPILER_IS_GNUCXX OR ${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
     SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
     SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
 
     set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -bind_at_load -m64")
     set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -bind_at_load -m64")
  endif()
endif()

message("compileflags: ${COMPILEFLAGS}")
add_definitions(${COMPILEFLAGS})


#dont change this destination. the GUI depends on it to guess the buildsystem

FILE(GLOB HEADERS "${SOURCE_DIR}/*.h*")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/)
list(APPEND deb_add_dirs "include/amdis")

set(include_sub_dirs "config" "parallel" "nonlin" "est" "expressions" 
					 "matrix_vector" "matrix_vector/expr" "operations" 
					 "traits" "utility" "time" "solver" "solver/details"
					 "solver/itl" "io" "io/detail")
				 
FOREACH(sub_dir ${include_sub_dirs})
	FILE(GLOB HEADERS "${SOURCE_DIR}/${sub_dir}/*.h*")
	INSTALL(FILES ${HEADERS} DESTINATION include/amdis/${sub_dir}/)
	list(APPEND deb_add_dirs "include/amdis/${sub_dir}")
ENDFOREACH(sub_dir)

FILE(GLOB HEADERS "${MUPARSER_SOURCE_DIR}/include/*.h")
INSTALL(FILES ${HEADERS}
	DESTINATION include/amdis/muparser)
list(APPEND deb_add_dirs "include/amdis/muparser")

list(APPEND deb_add_dirs "lib/amdis")
install(TARGETS amdis muparser
	DESTINATION lib/amdis/ )

configure_file(${AMDIS_SOURCE_DIR}/cmake/AMDISConfig.cmake.in
		${AMDIS_BINARY_DIR}/AMDISConfig.cmake
		@ONLY)

INSTALL(FILES ${AMDIS_BINARY_DIR}/AMDISConfig.cmake 
	DESTINATION share/amdis/)
list(APPEND deb_add_dirs "share/amdis/")

FILE(GLOB CMAKE_FILES "${AMDIS_SOURCE_DIR}/cmake/*.cmake")
INSTALL(FILES ${CMAKE_FILES}
	DESTINATION share/amdis/)

#install mtl
INSTALL(DIRECTORY ${AMDIS_SOURCE_DIR}/lib/mtl4/
	DESTINATION include/amdis/mtl4/
	FILES_MATCHING PATTERN "*.hpp"
	PATTERN ".svn" EXCLUDE
	PATTERN ".svn/*" EXCLUDE
	PATTERN "mtl4/libs" EXCLUDE
	PATTERN "mtl4/extern" EXCLUDE
       )
FILE(GLOB_RECURSE mtl4_files 
		RELATIVE ${AMDIS_SOURCE_DIR}/lib/mtl4/boost/ 
		${AMDIS_SOURCE_DIR}/lib/mtl4/boost/*.hpp)
FOREACH(mtl4_file in ${mtl4_files})
	get_filename_component(mtl4File_path ${mtl4_file} PATH)
	list(APPEND deb_add_dirs "include/amdis/mtl4/boost/${mtl4File_path}")
ENDFOREACH(mtl4_file)

#install documentation
FILE(GLOB DOCS "${AMDIS_SOURCE_DIR}/doc/*.pdf")
INSTALL(FILES  ${DOCS} DESTINATION share/amdis/doc)
list(APPEND deb_add_dirs "share/amdis/doc")

list(REMOVE_DUPLICATES deb_add_dirs)

ENABLE_TESTING()
INCLUDE(CTest)
INCLUDE(Dart)
if(NOT CTEST_MODEL)
        set(CTEST_MODEL "Experimental")
endif()
if( NOT N_TEST_PROC) 
	set(N_TEST_PROC "8")
endif()
configure_file(cmake/myMake.sh.in myMake.sh @ONLY)
configure_file(cmake/CTestScript.cmake.in CTestScript.cmake @ONLY)
