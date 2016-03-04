
if (ENABLE_EXTENSIONS)
    option(ENABLE_BASE_PROBLEMS "Use base_problems" true)
    mark_as_advanced(ENABLE_BASE_PROBLEMS)

    find_path(EXTENSIONS_DIR NAMES Helpers.h
        HINTS ${BASE_DIR}/../extensions
        DOC "Path to AMDiS extensions.")
    if (EXTENSIONS_DIR)
        if (NOT EXISTS ${EXTENSIONS_DIR}/Helpers.h OR NOT EXISTS ${EXTENSIONS_DIR}/ExtendedProblemStat.h)
            message(FATAL_ERROR "Wrong extensions directory! Directory must contain the files 'Helpers.h' and 'ExtendedProblemStat.h'")
        endif ()

        add_library(amdis_extensions
            ${EXTENSIONS_DIR}/Helpers.cc
            ${EXTENSIONS_DIR}/BackgroundMesh.cc
            ${EXTENSIONS_DIR}/GeometryTools.cc
            ${EXTENSIONS_DIR}/POperators.cc
            ${EXTENSIONS_DIR}/SingularDirichletBC2.cc
            ${EXTENSIONS_DIR}/time/ExtendedRosenbrockStationary.cc
            ${EXTENSIONS_DIR}/pugixml/src/pugixml.cpp
        )
        add_library(AMDiS::extensions ALIAS amdis_extensions)

        target_compile_definitions(amdis_extensions PRIVATE HAVE_EXTENSIONS=1)
        target_include_directories(amdis_extensions PRIVATE
            ${EXTENSIONS_DIR}
            ${EXTENSIONS_DIR}/time
            ${EXTENSIONS_DIR}/nanoflann
            ${EXTENSIONS_DIR}/pugixml/src
        )

        set(INSTALL_SUBDIRS . time preconditioner nanoflann)

        if (ENABLE_SEQ_PETSC)
            target_sources(amdis_extensions PRIVATE
                ${EXTENSIONS_DIR}/preconditioner/PetscPreconPfc.cc
                ${EXTENSIONS_DIR}/preconditioner/PetscPreconPfcDiag.cc
                ${EXTENSIONS_DIR}/preconditioner/PetscPreconCahnHilliard.cc)
        endif (ENABLE_SEQ_PETSC)

        if (ENABLE_PARALLEL_DOMAIN)
            target_sources(amdis_extensions PRIVATE
                ${EXTENSIONS_DIR}/preconditioner/PetscSolverPfc.cc
                ${EXTENSIONS_DIR}/preconditioner/PetscSolverPfc_diag.cc)
            target_link_libraries(amdis_extensions amdis_parallel)
        endif (ENABLE_PARALLEL_DOMAIN)

      if (ENABLE_BASE_PROBLEMS)
            if (ENABLE_REINIT)
                target_sources(amdis_extensions PRIVATE
                    ${EXTENSIONS_DIR}/base_problems/CahnHilliard.cc
                    ${EXTENSIONS_DIR}/base_problems/CahnHilliard_RB.cc
                    ${EXTENSIONS_DIR}/base_problems/CahnHilliardNavierStokes.cc
                    ${EXTENSIONS_DIR}/base_problems/CahnHilliardNavierStokes_RB.cc
                    ${EXTENSIONS_DIR}/base_problems/CahnHilliardNavierStokes_TwoPhase.cc
                    ${EXTENSIONS_DIR}/base_problems/CahnHilliardNavierStokes_TwoPhase_RB.cc)
            endif (ENABLE_REINIT)
            target_sources(amdis_extensions PRIVATE
                ${EXTENSIONS_DIR}/base_problems/DiffuseDomainFsi.cc
                ${EXTENSIONS_DIR}/base_problems/LinearElasticity.cc
                ${EXTENSIONS_DIR}/base_problems/LinearElasticityPhase.cc
                ${EXTENSIONS_DIR}/base_problems/NavierStokesCahnHilliard.cc
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
                ${EXTENSIONS_DIR}/base_problems/QuasiCrystal_RB.cc)
                # 		  ${EXTENSIONS_DIR}/base_problems/NavierStokes_Chorin.cc
                # 		  ${EXTENSIONS_DIR}/base_problems/NavierStokesPhase_Chorin.cc
                # 		  ${EXTENSIONS_DIR}/base_problems/VacancyPhaseFieldCrystal.cc

            target_compile_definitions(amdis_extensions PRIVATE HAVE_BASE_PROBLEMS=1)
            target_include_directories(amdis_extensions PRIVATE ${EXTENSIONS_DIR}/base_problems)

            list(APPEND INSTALL_SUBDIRS base_problems)
        endif (ENABLE_BASE_PROBLEMS)

        target_link_libraries(amdis_extensions amdis amdis_reinit)

        # specify how to install this target:
        # -----------------------------------

        foreach (SUBDIR ${INSTALL_SUBDIRS})
            file(GLOB HEADERS "${EXTENSIONS_DIR}/${SUBDIR}/*.h*")
            install(FILES ${HEADERS} DESTINATION include/amdis/extensions/${SUBDIR}/)
        endforeach ()

        file(GLOB HEADERS "${EXTENSIONS_DIR}/pugixml/src/*.hpp")
        install(FILES ${HEADERS} DESTINATION include/amdis/extensions/pugixml/)

        install(TARGETS amdis_extensions DESTINATION lib/amdis/)

    endif (EXTENSIONS_DIR)
endif (ENABLE_EXTENSIONS)
