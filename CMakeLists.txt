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
if(NOT CurrentRevision)
	find_package(Subversion)
	if(Subversion_FOUND)
	  Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} svn_info)
	  #message("svn-rev: ${svn_info_WC_REVISION}")
	  set(CurrentRevision "0.9.${svn_info_WC_REVISION}")
	endif()
endif()
if(NOT AMDIS_MAJOR)
  set(AMDIS_MAJOR "0")
endif()
if(NOT AMDIS_MINOR)
  set(AMDIS_MINOR "9")
endif()
set(COMPILEFLAGS "-DAMDIS_VERSION=${CurrentRevision}")

#define the build type, empty can be everything and nothing
if(CMAKE_BUILD_TYPE STREQUAL "")
	set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_CXX_COMPILER MATCHES ".*icpc")
	Message("Compiler: Intel")
	set(CMAKE_CXX_FLAGS "-diag-disable 654 -diag-disable 858") 
endif()


SET(ENABLE_PARALLEL_DOMAIN "OFF" CACHE STRING "use parallel domain decomposition. please set to one of: PMTL, PETSC, OFF" )
option(USE_PETSC_DEV false)
option(ENABLE_ZOLTAN false)
option(ENABLE_SEQ_PETSC "Use sequential PETSc solvers" false)
option(ENABLE_UMFPACK "Use of UMFPACK solver" false)
option(ENABLE_HYPRE "Use HYPRE AMG solver" false)
option(ENABLE_PNG "use png reader/writer" false)
option(ENABLE_BDDCML "Use of BDDCML library" false)
option(ENABLE_EXTENSIONS "Use extensions" false)
option(ENABLE_OPENMP "Use OpenMP" false)
option(ENABLE_COMPRESSION "Use output compression" false)
option(ENABLE_OUTPUT "AMDiS output printing, disable only for debugging!" true)


mark_as_advanced(ENABLE_OUTPUT)
mark_as_advanced(ENABLE_BDDCML)
mark_as_advanced(USE_PETSC_DEV)


find_package(Boost 1.44 REQUIRED)
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
	      ${SOURCE_DIR}/Bubble.cc
	      ${SOURCE_DIR}/Cholesky.cc
	      ${SOURCE_DIR}/CoarseningManager.cc
	      ${SOURCE_DIR}/CoarseningManager1d.cc
	      ${SOURCE_DIR}/CoarseningManager2d.cc
	      ${SOURCE_DIR}/CoarseningManager3d.cc
	      ${SOURCE_DIR}/ComponentTraverseInfo.cc
	      ${SOURCE_DIR}/CouplingIterationInterface.cc
	      ${SOURCE_DIR}/CreatorMap.cc
	      ${SOURCE_DIR}/DOFAdmin.cc
	      ${SOURCE_DIR}/DOFIndexed.cc
	      ${SOURCE_DIR}/DOFMatrix.cc
	      ${SOURCE_DIR}/DOFVector.cc
	      ${SOURCE_DIR}/Debug.cc
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
	      ${SOURCE_DIR}/GlobalDOFNumbering.cc
	      ${SOURCE_DIR}/GlobalElementNumbering.cc
	      ${SOURCE_DIR}/Initfile.cc
	      ${SOURCE_DIR}/Lagrange.cc
	      ${SOURCE_DIR}/LeafData.cc
	      ${SOURCE_DIR}/Line.cc
	      ${SOURCE_DIR}/MacroElement.cc
	      ${SOURCE_DIR}/Marker.cc
	      ${SOURCE_DIR}/MatrixVector.cc
# 	      ${SOURCE_DIR}/Test_MatrixVectorOperations.cc
	      ${SOURCE_DIR}/Mesh.cc
	      ${SOURCE_DIR}/MeshStructure.cc
	      ${SOURCE_DIR}/Operator.cc
	      ${SOURCE_DIR}/OperatorTerm.cc
	      ${SOURCE_DIR}/Parametric.cc
	      ${SOURCE_DIR}/PeriodicBC.cc
	      #${SOURCE_DIR}/ProblemImplicit.cc
	      ${SOURCE_DIR}/ProblemInstat.cc
	      ${SOURCE_DIR}/ProblemInterpol.cc
	      ${SOURCE_DIR}/ProblemStat.cc
	      ${SOURCE_DIR}/ProblemStatDbg.cc
	      ${SOURCE_DIR}/Projection.cc
	      ${SOURCE_DIR}/QPInfo.cc
	      ${SOURCE_DIR}/QPsiPhi.cc
	      ${SOURCE_DIR}/Quadrature.cc
	      ${SOURCE_DIR}/RCNeighbourList.cc
	      ${SOURCE_DIR}/Recovery.cc
	      ${SOURCE_DIR}/RefinementManager.cc
	      ${SOURCE_DIR}/RefinementManager1d.cc
	      ${SOURCE_DIR}/RefinementManager2d.cc
	      ${SOURCE_DIR}/RefinementManager3d.cc
	      ${SOURCE_DIR}/RobinBC.cc
	      ${SOURCE_DIR}/ScalableQuadrature.cc
	      ${SOURCE_DIR}/SecondOrderAssembler.cc
	      ${SOURCE_DIR}/SecondOrderTerm.cc
	      ${SOURCE_DIR}/Serializer.cc
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
	      ${SOURCE_DIR}/ZeroOrderTerm.cc
	      ${SOURCE_DIR}/est/Estimator.cc
	      ${SOURCE_DIR}/est/RecoveryEstimator.cc
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
	      ${SOURCE_DIR}/io/VtkWriter.cc
	      ${SOURCE_DIR}/io/VtkVectorWriter.cc
              ${SOURCE_DIR}/io/detail/VtkWriter.cc
	      ${SOURCE_DIR}/nonlin/ProblemNonLin.cc
	      ${SOURCE_DIR}/solver/SolverMatrix.cc
	      ${SOURCE_DIR}/time/RosenbrockAdaptInstationary.cc
	      ${SOURCE_DIR}/time/RosenbrockMethod.cc
	      ${SOURCE_DIR}/time/RosenbrockStationary.cc
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

if(ENABLE_COMPRESSION)
  list(APPEND COMPILEFLAGS "-DHAVE_COMPRESSION=1")
endif(ENABLE_COMPRESSION)

if(NOT ENABLE_OUTPUT)
       message(WARNING "AMDiS cout output disabled!")
       list(APPEND COMPILEFLAGS "-DSUPPRESS_OUTPUT")
endif(NOT ENABLE_OUTPUT)


if(ENABLE_BDDCML)
	SET(BDDCML_LINK_LIST "" CACHE STRING "Further libraries to be linked with BDDCML")
		
	find_file(BDDCML_H bddcml_interface_c.h
			   HINTS ENV CPATH
			   DOC "Header bddcml_interface_c.h for the BDDCML library.")

	if(BDDCML_H)
		get_filename_component(BDDCML_PATH ${BDDCML_H} PATH)
		include_directories(${BDDCML_PATH})		
		list(APPEND COMPILEFLAGS "-DHAVE_BDDC_ML=1")
		list(APPEND COMPILEFLAGS "-DAdd_")
	else()
		message(FATAL_ERROR "Could not find BDDCML headers.")
	endif()

	find_library(BDDCML_LIB bddcml
				HINTS ENV LIBRARY_PATH
				DOC "BDDCML library")

	if(BDDCML_LIB)
		list(APPEND AMDIS_LIBS ${BDDCML_LIB})
	else()
		message(FATAL_ERROR "Could not find the BDDCML library")
	endif()


	find_file(BLOPEX_PATH blopex
			      HINTS ENV CPATH
			      DOC "Path to BLOPEX library.")
endif(ENABLE_BDDCML)


if(ENABLE_EXTENSIONS)
  option(ENABLE_BASE_PROBLEMS "Use base_problems" false)
  find_path(EXTENSIONS_DIR NAMES ${SOURCE_DIR}/../../extensions/Helpers.h DOC "Path to AMDiS extensions.")
  if (EXTENSIONS_DIR)
    if (NOT EXISTS ${EXTENSIONS_DIR}/Helpers.h OR NOT EXISTS ${EXTENSIONS_DIR}/base_problems/BaseProblem.h)
      message(FATAL_ERROR "Wrong extensions directory! Directory must contain the 'Helpers.h' and the subdirectory 'base_problems'")
    else()
  
      SET(EXTENSIONS_SRC
	      ${EXTENSIONS_DIR}/Helpers.cc
	      ${EXTENSIONS_DIR}/BackgroundMesh.cc
	      ${EXTENSIONS_DIR}/GeometryTools.cc
	      ${EXTENSIONS_DIR}/POperators.cc
	      ${EXTENSIONS_DIR}/SingularDirichletBC2.cc
	      ${EXTENSIONS_DIR}/time/ExtendedRosenbrockStationary.cc
	      ${EXTENSIONS_DIR}/pugixml/src/pugixml.cpp)
	      
      if(ENABLE_SEQ_PETSC)
	list(APPEND EXTENSIONS_SRC
	      ${EXTENSIONS_DIR}/preconditioner/PetscPreconPfc.cc
	      ${EXTENSIONS_DIR}/preconditioner/PetscPreconPfcDiag.cc
	      ${EXTENSIONS_DIR}/preconditioner/PetscPreconCahnHilliard.cc)
      endif(ENABLE_SEQ_PETSC)
      
      if(ENABLE_PARALLEL_DOMAIN)
	list(APPEND EXTENSIONS_SRC
# 	      ${EXTENSIONS_DIR}/preconditioner/PetscSolverNavierStokes2.cc
	      ${EXTENSIONS_DIR}/preconditioner/PetscSolverPfc.cc
	      ${EXTENSIONS_DIR}/preconditioner/PetscSolverPfc_diag.cc)
      endif(ENABLE_PARALLEL_DOMAIN)
      
      list(APPEND COMPILEFLAGS "-DHAVE_EXTENSIONS=1")
      list(APPEND AMDIS_INCLUDE_DIRS ${EXTENSIONS_DIR})
      list(APPEND AMDIS_INCLUDE_DIRS ${EXTENSIONS_DIR}/time)
      list(APPEND AMDIS_INCLUDE_DIRS ${EXTENSIONS_DIR}/nanoflann)
      list(APPEND AMDIS_INCLUDE_DIRS ${EXTENSIONS_DIR}/pugixml/src)
      if(WIN32)
	add_library(extensions STATIC ${EXTENSIONS_SRC})
      else()
	add_library(extensions SHARED ${EXTENSIONS_SRC})
      endif()
      list(APPEND AMDIS_LIBS extensions)

      FILE(GLOB HEADERS "${EXTENSIONS_DIR}/*.h")
      INSTALL(FILES ${HEADERS}
	      DESTINATION include/amdis/extensions/)
      FILE(GLOB HEADERS "${EXTENSIONS_DIR}/*.hh")
      INSTALL(FILES ${HEADERS}
	      DESTINATION include/amdis/extensions/)

      list(APPEND deb_add_dirs "include/amdis/extensions")

      FILE(GLOB HEADERS "${EXTENSIONS_DIR}/time/*.h")
      INSTALL(FILES ${HEADERS}
	      DESTINATION include/amdis/extensions/time/)
      FILE(GLOB HEADERS "${EXTENSIONS_DIR}/time/*.hh")
      INSTALL(FILES ${HEADERS}
	      DESTINATION include/amdis/extensions/time/)
	      
      FILE(GLOB HEADERS "${EXTENSIONS_DIR}/preconditioner/*.h")
      INSTALL(FILES ${HEADERS}
	      DESTINATION include/amdis/extensions/preconditioner/)
      FILE(GLOB HEADERS "${EXTENSIONS_DIR}/preconditioner/*.hh")
      INSTALL(FILES ${HEADERS}
	      DESTINATION include/amdis/extensions/preconditioner/)

      list(APPEND deb_add_dirs "include/amdis/extensions/time")

      FILE(GLOB HEADERS "${EXTENSIONS_DIR}/nanoflann/*.hpp")
      INSTALL(FILES ${HEADERS}
	      DESTINATION include/amdis/extensions/nanoflann/)

      list(APPEND deb_add_dirs "include/amdis/extensions/nanoflann")

      FILE(GLOB HEADERS "${EXTENSIONS_DIR}/pugixml/src/*.hpp")
      INSTALL(FILES ${HEADERS}
	      DESTINATION include/amdis/extensions/pugixml/)

      list(APPEND deb_add_dirs "include/amdis/extensions/pugixml")

      install(TARGETS extensions
	      DESTINATION lib/amdis/ )

# === baseProblems ====================================================================
      
      if(ENABLE_BASE_PROBLEMS)
	  SET(BASE_PROBLEMS_SRC
		  ${EXTENSIONS_DIR}/base_problems/CahnHilliard.cc
		  ${EXTENSIONS_DIR}/base_problems/CahnHilliard_RB.cc
		  ${EXTENSIONS_DIR}/base_problems/CahnHilliardNavierStokes.cc
		  ${EXTENSIONS_DIR}/base_problems/CahnHilliardNavierStokes_RB.cc
		  ${EXTENSIONS_DIR}/base_problems/CahnHilliardNavierStokes_TwoPhase.cc
		  ${EXTENSIONS_DIR}/base_problems/CahnHilliardNavierStokes_TwoPhase_RB.cc
		  ${EXTENSIONS_DIR}/base_problems/DiffuseDomainFsi.cc
		  ${EXTENSIONS_DIR}/base_problems/LinearElasticity.cc
		  ${EXTENSIONS_DIR}/base_problems/LinearElasticityPhase.cc
# 		  ${EXTENSIONS_DIR}/base_problems/NavierStokes_Chorin.cc
		  ${EXTENSIONS_DIR}/base_problems/NavierStokesCahnHilliard.cc
# 		  ${EXTENSIONS_DIR}/base_problems/NavierStokesPhase_Chorin.cc
		  ${EXTENSIONS_DIR}/base_problems/NavierStokesPhase_TaylorHood.cc
		  ${EXTENSIONS_DIR}/base_problems/NavierStokes_TaylorHood.cc
		  ${EXTENSIONS_DIR}/base_problems/NavierStokes_TaylorHood_RB.cc
		  ${EXTENSIONS_DIR}/base_problems/NavierStokes_TH_MultiPhase.cc
		  ${EXTENSIONS_DIR}/base_problems/NavierStokes_TH_MultiPhase_RB.cc
		  ${EXTENSIONS_DIR}/base_problems/PhaseFieldCrystal.cc
		  ${EXTENSIONS_DIR}/base_problems/PhaseFieldCrystal_Phase.cc
		  ${EXTENSIONS_DIR}/base_problems/PhaseFieldCrystal_RB.cc
		  ${EXTENSIONS_DIR}/base_problems/PolarizationField.cc
		  ${EXTENSIONS_DIR}/base_problems/QuasiCrystal.cc
		  ${EXTENSIONS_DIR}/base_problems/QuasiCrystal_RB.cc
# 		  ${EXTENSIONS_DIR}/base_problems/VacancyPhaseFieldCrystal.cc
		  )
	list(APPEND COMPILEFLAGS "-DHAVE_BASE_PROBLEMS=1")
	list(APPEND AMDIS_INCLUDE_DIRS ${EXTENSIONS_DIR}/base_problems)
	if(WIN32)
	  add_library(base_problems STATIC ${BASE_PROBLEMS_SRC})
	else()
	  add_library(base_problems SHARED ${BASE_PROBLEMS_SRC})
	endif()
	list(APPEND AMDIS_LIBS base_problems)

	FILE(GLOB HEADERS "${EXTENSIONS_DIR}/base_problems/*.h")
	INSTALL(FILES ${HEADERS}
		DESTINATION include/amdis/base_problems/)
	FILE(GLOB HEADERS "${EXTENSIONS_DIR}/base_problems/*.hh")
	INSTALL(FILES ${HEADERS}
		DESTINATION include/amdis/base_problems/)

	list(APPEND deb_add_dirs "include/amdis/base_problems")

	install(TARGETS base_problems
		DESTINATION lib/amdis/ )
    
      endif(ENABLE_BASE_PROBLEMS)
    endif()
  endif()
endif(ENABLE_EXTENSIONS)


if(ENABLE_OPENMP)
  if(CMAKE_CXX_COMPILER MATCHES ".*icpc")
	list(APPEND COMPILEFLAGS "-openmp")
  else()
        list(APPEND COMPILEFLAGS "-fopenmp")
  endif()
endif(ENABLE_OPENMP)

SET(COMPOSITE_SOURCE_DIR ${SOURCE_DIR}/compositeFEM)
SET(COMPOSITE_FEM_SRC ${COMPOSITE_SOURCE_DIR}/CFE_Integration.cc 
		      ${COMPOSITE_SOURCE_DIR}/CFE_NormAndErrorFcts.cc 
		      ${COMPOSITE_SOURCE_DIR}/CompositeFEMMethods.cc 
		      ${COMPOSITE_SOURCE_DIR}/CompositeFEMOperator.cc 
		      ${COMPOSITE_SOURCE_DIR}/LevelSetAdaptMesh.cc 
		      ${COMPOSITE_SOURCE_DIR}/PenaltyOperator.cc 
		      ${COMPOSITE_SOURCE_DIR}/ElementLevelSet.cc 
		      ${COMPOSITE_SOURCE_DIR}/SubPolytope.cc 
		      ${COMPOSITE_SOURCE_DIR}/SubElementAssembler.cc )

set(REINIT_SOURCE_DIR ${SOURCE_DIR}/reinit)
file(GLOB REINIT_SRC ${REINIT_SOURCE_DIR}/*.cc)
include_directories(${REINIT_SOURCE_DIR})


#muparser includes
set(MUPARSER_SOURCE_DIR ${AMDIS_SOURCE_DIR}/lib/muparser_v134)
file(GLOB MUPARSER_SRC ${MUPARSER_SOURCE_DIR}/src/*.cpp)
list(APPEND AMDIS_INCLUDE_DIRS ${MUPARSER_SOURCE_DIR}/include)

#mtl4 includes
list(APPEND AMDIS_INCLUDE_DIRS ${MTL_INCLUDE_DIR})
#include_directories(${MTL_INCLUDE_DIR})
list(APPEND AMDIS_INCLUDE_DIRS ${SOURCE_DIR})

#include_directories(${SOURCE_DIR})
include_directories(${AMDIS_INCLUDE_DIRS})

if(WIN32)
  add_library(amdis STATIC ${AMDIS_SRC} ${PARALLEL_DOMAIN_AMDIS_SRC})
  add_library(compositeFEM STATIC ${COMPOSITE_FEM_SRC})
  add_library(reinit STATIC ${REINIT_SRC})
  add_library(muparser STATIC ${MUPARSER_SRC})
else()
  add_library(amdis SHARED ${AMDIS_SRC} ${PARALLEL_DOMAIN_AMDIS_SRC})
  add_library(compositeFEM SHARED ${COMPOSITE_FEM_SRC})
  add_library(reinit SHARED ${REINIT_SRC})
  add_library(muparser SHARED ${MUPARSER_SRC})
endif()

#target_link_libraries(compositeFEM amdis)
#target_link_libraries(reinit amdis)
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
#SET_TARGET_PROPERTIES(amdis PROPERTIES COMPILE_FLAGS "${COMPILEFLAGS}")
add_definitions(${COMPILEFLAGS})

if(ENABLE_MARMOT)
	Message("please set marmotcc manually")
endif(ENABLE_MARMOT)

#dont change this destination. the GUI depends on it to guess the buildsystem

FILE(GLOB HEADERS "${SOURCE_DIR}/*.h")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/)
list(APPEND deb_add_dirs "include/amdis")

FILE(GLOB HEADERS "${SOURCE_DIR}/config/*.h*")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/config/)
list(APPEND deb_add_dirs "include/amdis/config")

FILE(GLOB HEADERS "${SOURCE_DIR}/*.hh")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/)

FILE(GLOB HEADERS "${SOURCE_DIR}/*.hpp")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/)

FILE(GLOB HEADERS "${SOURCE_DIR}/parallel/*.h")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/parallel/)
list(APPEND deb_add_dirs "include/amdis/parallel")

FILE(GLOB HEADERS "${SOURCE_DIR}/nonlin/*.h")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/nonlin/)
list(APPEND deb_add_dirs "include/amdis/nonlin")

FILE(GLOB HEADERS "${SOURCE_DIR}/est/*.h")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/est/)
list(APPEND deb_add_dirs "include/amdis/est")

FILE(GLOB HEADERS "${SOURCE_DIR}/expressions/*.h*")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/expressions/)
list(APPEND deb_add_dirs "include/amdis/expressions")

FILE(GLOB HEADERS "${SOURCE_DIR}/operations/*.h*")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/operations/)
list(APPEND deb_add_dirs "include/amdis/operations")

FILE(GLOB HEADERS "${SOURCE_DIR}/traits/*.h*")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/traits/)
list(APPEND deb_add_dirs "include/amdis/traits")

FILE(GLOB HEADERS "${SOURCE_DIR}/utility/*.h*")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/utility/)
list(APPEND deb_add_dirs "include/amdis/utility")

FILE(GLOB HEADERS "${SOURCE_DIR}/time/*.h")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/time/)
list(APPEND deb_add_dirs "include/amdis/time")

FILE(GLOB HEADERS "${SOURCE_DIR}/solver/*.h")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/solver)
list(APPEND deb_add_dirs "include/amdis/solver")

FILE(GLOB HEADERS "${SOURCE_DIR}/solver/*.hh")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/solver)
list(APPEND deb_add_dirs "include/amdis/solver")

FILE(GLOB HEADERS "${SOURCE_DIR}/solver/details/*.h")
INSTALL(FILES ${HEADERS} 
        DESTINATION include/amdis/solver/details)
list(APPEND deb_add_dirs "include/amdis/solver/details")

FILE(GLOB HEADERS "${SOURCE_DIR}/solver/itl/*.hpp")
INSTALL(FILES ${HEADERS} 
	DESTINATION include/amdis/solver/itl)
list(APPEND deb_add_dirs "include/amdis/solver/itl")

FILE(GLOB HEADERS "${SOURCE_DIR}/io/*.h")
INSTALL(FILES ${HEADERS}
	DESTINATION include/amdis/io)
list(APPEND deb_add_dirs "include/amdis/io")

FILE(GLOB HEADERS "${SOURCE_DIR}/io/*.hh")
INSTALL(FILES ${HEADERS}
	DESTINATION include/amdis/io)
	
FILE(GLOB HEADERS "${SOURCE_DIR}/io/detail/*h")
INSTALL(FILES ${HEADERS}
        DESTINATION include/amdis/io/detail)
list(APPEND deb_add_dirs "include/amdis/io/detail")        

FILE(GLOB HEADERS "${REINIT_SOURCE_DIR}/*.h")
INSTALL(FILES ${HEADERS}
	DESTINATION include/amdis/reinit)
list(APPEND deb_add_dirs "include/amdis/reinit")

FILE(GLOB HEADERS "${COMPOSITE_SOURCE_DIR}/*.h")
INSTALL(FILES ${HEADERS}
	DESTINATION include/amdis/compositeFEM)
list(APPEND deb_add_dirs "include/amdis/compositeFEM")

FILE(GLOB HEADERS "${MUPARSER_SOURCE_DIR}/include/*.h")
INSTALL(FILES ${HEADERS}
	DESTINATION include/amdis/muparser)
list(APPEND deb_add_dirs "include/amdis/muparser")

list(APPEND deb_add_dirs "lib/amdis")
install(TARGETS amdis compositeFEM reinit muparser
	DESTINATION lib/amdis/ )

configure_file(${AMDIS_SOURCE_DIR}/AMDISConfig.cmake.in
		${AMDIS_BINARY_DIR}/AMDISConfig.cmake
		@ONLY
	      )

INSTALL(FILES ${AMDIS_BINARY_DIR}/AMDISConfig.cmake 
	DESTINATION share/amdis/)
list(APPEND deb_add_dirs "share/amdis/")

# INSTALL(FILES ${AMDIS_SOURCE_DIR}/AMDISUse.cmake 
# 	DESTINATION share/amdis/)
FILE(GLOB CMAKE_FILES "${AMDIS_SOURCE_DIR}/*.cmake")
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
foreach(mtl4_file in ${mtl4_files})
 get_filename_component(mtl4File_path ${mtl4_file} PATH)
 list(APPEND deb_add_dirs "include/amdis/mtl4/boost/${mtl4File_path}")
endforeach(mtl4_file)

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
configure_file(myMake.sh.in myMake.sh @ONLY)
configure_file(CTestScript.cmake.in CTestScript.cmake @ONLY)


# PACKAGES
# ========================================================
set(CPACK_PACKAGE_NAME "AMDIS")
set(CPACK_PACKAGE_CONTACT "Andreas.Naumann@tu-dresden.de")
set(CPACK_PACKAGE_VERSION "${CurrentRevision}")
set(CPACK_GENERATOR "RPM" "DEB" "TGZ")
set(CPACK_INSTALL_PREFIX "/usr/")
set(CPACK_RESOURCE_FILE_LICENSE "${AMDIS_SOURCE_DIR}/license.opensource.txt")

# library depends on boost, but boost libraries have version-number in name on some systems.
# so we have to find a workaround to enforce allowed boost versions
set(boost_ver ${Boost_MINOR_VERSION})
math(EXPR boost_ver_start "${Boost_MINOR_VERSION}+1")
math(EXPR boost_ver_stop "${Boost_MINOR_VERSION}+10")

set(boost_depends "libboost${Boost_MAJOR_VERSION}.${boost_ver}-dev (>= ${Boost_MAJOR_VERSION}.${boost_ver})")
set(boost_iostreams_depends "libboost-iostreams${Boost_MAJOR_VERSION}.${boost_ver}-dev (>= ${Boost_MAJOR_VERSION}.${boost_ver})")
set(boost_system_depends "libboost-system${Boost_MAJOR_VERSION}.${boost_ver}-dev (>= ${Boost_MAJOR_VERSION}.${boost_ver})")
set(boost_filesystem_depends "libboost-filesystem${Boost_MAJOR_VERSION}.${boost_ver}-dev (>= ${Boost_MAJOR_VERSION}.${boost_ver})")

foreach(boost_ver RANGE ${boost_ver_start} ${boost_ver_stop})
  set(boost_depends "${boost_depends} | libboost${Boost_MAJOR_VERSION}.${boost_ver}-dev (>= ${Boost_MAJOR_VERSION}.${boost_ver})")
  set(boost_iostreams_depends "${boost_iostreams_depends} | libboost-iostreams${Boost_MAJOR_VERSION}.${boost_ver}-dev (>= ${Boost_MAJOR_VERSION}.${boost_ver})")
  set(boost_system_depends "${boost_system_depends} | libboost-system${Boost_MAJOR_VERSION}.${boost_ver}-dev (>= ${Boost_MAJOR_VERSION}.${boost_ver})")
  set(boost_filesystem_depends "${boost_filesystem_depends} | libboost-filesystem${Boost_MAJOR_VERSION}.${boost_ver}-dev (>= ${Boost_MAJOR_VERSION}.${boost_ver})")
endforeach(boost_ver)

configure_file(${AMDIS_SOURCE_DIR}/preinst.in 
  ${AMDIS_BINARY_DIR}/preinst 
  @ONLY)

configure_file(${AMDIS_SOURCE_DIR}/postrm.in 
  ${AMDIS_BINARY_DIR}/postrm 
  @ONLY)
  
#----Test if clang setup works----------------------------------------------------------------------
set(COMPILER_VERSION "")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  exec_program(${CMAKE_C_COMPILER} ARGS "-v" OUTPUT_VARIABLE _clang_version_info)
  string(REGEX REPLACE "^.*[ ]([0-9]+)\\.[0-9].*$" "\\1" CLANG_MAJOR "${_clang_version_info}")
  string(REGEX REPLACE "^.*[ ][0-9]+\\.([0-9]).*$" "\\1" CLANG_MINOR "${_clang_version_info}")
  set(COMPILER_VERSION "clang ${CLANG_MAJOR}.${CLANG_MINOR}")
endif()

#---Obtain the major and minor version of the GNU compiler-------------------------------------------
if (CMAKE_COMPILER_IS_GNUCXX)
  exec_program(${CMAKE_C_COMPILER} ARGS "-dumpversion" OUTPUT_VARIABLE _gcc_version_info)
  string(REGEX REPLACE "^([0-9]+).*$"                   "\\1" GCC_MAJOR ${_gcc_version_info})
  string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*$"          "\\1" GCC_MINOR ${_gcc_version_info})
  string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" GCC_PATCH ${_gcc_version_info})

  if(GCC_PATCH MATCHES "\\.+")
    set(GCC_PATCH "")
  endif()
  if(GCC_MINOR MATCHES "\\.+")
    set(GCC_MINOR "")
  endif()
  if(GCC_MAJOR MATCHES "\\.+")
    set(GCC_MAJOR "")
  endif()
  set(COMPILER_VERSION "gcc ${_gcc_version_info}")
endif()
  
set(CPACK_RPM_PACKAGE_REQUIRES "boost-devel >= ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}, ${RPM_DEPEND_STR}")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "${boost_depends},\\n ${boost_iostreams_depends},\\n ${boost_system_depends},\\n ${boost_filesystem_depends}")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "AMDiS (Adaptive Multi-Dimensional Simulations) - a C++ library to solve partial differential equations using adaptive finite elements. The main features of AMDiS are:
 . 
  - full adaptive meshes in 1D, 2D and 3D
  - adaptive time stepping strategies for time dependent equations
  - uses high performance Krylov subspace solvers from MTL4 and the direct solvers UMFPACK and Intel's Pardiso
  - parallelization of the overall finite element method with excelent scaling up to several thounds of processors
  - uses PETSc for solving parallel distributed linear systems
  - easy to use interface, no additional effort to run code in parallel!
 . 
  The library is compiled with ${COMPILER_VERSION} and compiled/linked with boost ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}.")

set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${AMDIS_BINARY_DIR}/preinst" "${AMDIS_BINARY_DIR}/postrm")

include(CPack)
