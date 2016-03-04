if (ENABLE_PARALLEL_DOMAIN)
    option(ENABLE_ZOLTAN "Add support for the Parallel Partitioning suite Zoltan" false)
    option(ENABLE_PARALLEL_SOLVERS "Add some problem dependent solver, e.g. Feti, Navier-Stokes and Cahn-Hilliard" true)
    mark_as_advanced(ENABLE_PARALLEL_SOLVERS)
    
    add_library(amdis_parallel INTERFACE)
    target_sources(amdis PRIVATE
        ${SOURCE_DIR}/parallel/DofComm.cc
        ${SOURCE_DIR}/parallel/CheckerPartitioner.cc
        ${SOURCE_DIR}/parallel/ElementObjectDatabase.cc
        ${SOURCE_DIR}/parallel/InteriorBoundary.cc
        ${SOURCE_DIR}/parallel/MeshDistributor.cc
        ${SOURCE_DIR}/parallel/MeshLevelData.cc
        ${SOURCE_DIR}/parallel/MeshManipulation.cc
        ${SOURCE_DIR}/parallel/MeshPartitioner.cc
        ${SOURCE_DIR}/parallel/MpiHelper.cc
        ${SOURCE_DIR}/parallel/ParallelDofMapping.cc
        ${SOURCE_DIR}/parallel/ParallelProblemStat.cc
        ${SOURCE_DIR}/parallel/ParallelSolver.cc
        ${SOURCE_DIR}/parallel/PeriodicMap.cc
        ${SOURCE_DIR}/parallel/ParMetisPartitioner.cc
        ${SOURCE_DIR}/parallel/StdMpi.cc
    )
    
    target_sources(amdis_debug INTERFACE
        ${SOURCE_DIR}/parallel/ParallelDebug.cc
    )
  
    target_compile_definitions(amdis_parallel INTERFACE
	   HAVE_PARALLEL_DOMAIN_AMDIS=1)
  
    # MPI is required
    find_package(MPI REQUIRED)
    if (MPI_FOUND)
        target_include_directories(amdis_parallel INTERFACE
            ${MPI_INCLUDE_PATH})
        target_compile_options(amdis_parallel INTERFACE
            ${MPI_COMPILE_FLAGS})
    endif (MPI_FOUND)
  
    # PETSc library is required
    set(PETSC_EXECUTABLE_RUNS ON)
    include(find_petsc)
    if (PETSc_FOUND)
        target_include_directories(amdis_parallel INTERFACE
            ${PETSC_DIR}/include 
            ${PETSC_DIR}/${PETSC_ARCH}/include)
        
        # parmetis is required
        find_file(PARMETIS_HEADER_FILE "parmetis.h" HINTS ${PETSC_DIR}/include )
        if (PARMETIS_HEADER_FILE)
            get_filename_component(PARMETIS_INCLUDE_PATH "${PARMETIS_HEADER_FILE}" PATH CACHE)
            target_include_directories(amdis_parallel INTERFACE ${PARMETIS_INCLUDE_PATH})
        else()
            message(FATAL_ERROR "Could not find ParMetis header file 'parmetis.h'!")
        endif (PARMETIS_HEADER_FILE)
        
        # add support for the zoltan library
        if (ENABLE_ZOLTAN)
            find_file(ZOLTAN_HEADER_FILE "zoltan_cpp.h" HINTS ${PETSC_DIR}/include)
            if (ZOLTAN_HEADER_FILE)
                get_filename_component(ZOLTAN_HEADER_DIR "${ZOLTAN_HEADER_FILE}" PATH CACHE)
                target_include_directories(amdis_parallel PRIINTERFACEVATE ${ZOLTAN_HEADER_DIR})
            else()
                message(FATAL_ERROR "Could not find Zoltan include file 'zoltan_cpp.h'!")
            endif(ZOLTAN_HEADER_FILE)
            
            target_compile_definitions(amdis_parallel INTERFACE HAVE_ZOLTAN=1)
            target_sources(amdis PRIVATE
                ${SOURCE_DIR}/parallel/ZoltanPartitioner.cc)
        endif (ENABLE_ZOLTAN)
        
        if (ENABLE_BDDCML)
            target_compile_definitions(amdis_parallel INTERFACE HAVE_BDDCML=1)
            target_sources(amdis PRIVATE
                ${SOURCE_DIR}/parallel/BddcMlSolver.cc)
        endif (ENABLE_BDDCML)
        
        # add some more source-files that need petsc
        target_sources(amdis PRIVATE
            ${SOURCE_DIR}/parallel/MatrixNnzStructure.cc
            ${SOURCE_DIR}/parallel/ParallelCoarseSpaceSolver.cc
            ${SOURCE_DIR}/parallel/PetscHelper.cc
            ${SOURCE_DIR}/parallel/PetscSolver.cc
            ${SOURCE_DIR}/parallel/PetscSolverGlobalMatrix.cc
            ${SOURCE_DIR}/parallel/PetscSolverGlobalBlockMatrix.cc
            ${SOURCE_DIR}/solver/PetscTypes.cc
        )
            
        if (ENABLE_PARALLEL_SOLVERS)
            target_sources(amdis PRIVATE
                ${SOURCE_DIR}/parallel/PetscSolverFeti.cc
                ${SOURCE_DIR}/parallel/PetscSolverFetiDebug.cc
                ${SOURCE_DIR}/parallel/PetscSolverFetiMonitor.cc
                ${SOURCE_DIR}/parallel/PetscSolverFetiOperators.cc
                ${SOURCE_DIR}/parallel/PetscSolverFetiTimings.cc
                ${SOURCE_DIR}/parallel/PetscSolverNavierStokes.cc
                ${SOURCE_DIR}/parallel/PetscSolverNSCH.cc
                ${SOURCE_DIR}/parallel/PetscSolverCahnHilliard2.cc
                ${SOURCE_DIR}/parallel/PetscSolverCahnHilliard.cc
                ${SOURCE_DIR}/parallel/PetscSolverSchur.cc
            )
            target_compile_definitions(amdis_parallel INTERFACE
                HAVE_PARALLEL_SOLVERS=1)
        endif (ENABLE_PARALLEL_SOLVERS)
            
        target_compile_definitions(amdis_parallel INTERFACE
            HAVE_PARALLEL_PETSC=1
            PETSC_VERSION=${PETSC_VERSION})
            
        target_link_libraries(amdis amdis_parallel ${PETSC_LIBRARIES} blas lapack)
    endif (PETSc_FOUND)
    
    # specify how to install this target:
    # -----------------------------------
    
    file(GLOB AMDIS_PARALLEL_HEADERS "${SOURCE_DIR}/parallel/*.h")
    install(FILES ${AMDIS_PARALLEL_HEADERS} DESTINATION include/amdis/parallel/)
    
    install(FILES 
        ${BASE_DIR}/cmake3/ResolveCompilerPaths.cmake 
        ${BASE_DIR}/cmake3/FindPackageMultipass.cmake 
        ${BASE_DIR}/cmake3/find_petsc.cmake
        DESTINATION share/amdis/)
endif (ENABLE_PARALLEL_DOMAIN)
