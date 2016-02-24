# specify the target and requirements for the composite fem library

SET(COMPOSITE_SOURCE_DIR ${SOURCE_DIR}/compositeFEM)
add_library(amdis_compositeFEM
    ${COMPOSITE_SOURCE_DIR}/CFE_Integration.cc
    ${COMPOSITE_SOURCE_DIR}/CFE_NormAndErrorFcts.cc
    ${COMPOSITE_SOURCE_DIR}/CompositeFEMMethods.cc
    ${COMPOSITE_SOURCE_DIR}/CompositeFEMOperator.cc
    ${COMPOSITE_SOURCE_DIR}/LevelSetAdaptMesh.cc
    ${COMPOSITE_SOURCE_DIR}/PenaltyOperator.cc
    ${COMPOSITE_SOURCE_DIR}/ElementLevelSet.cc
    ${COMPOSITE_SOURCE_DIR}/SubPolytope.cc
    ${COMPOSITE_SOURCE_DIR}/SubElementAssembler.cc)
add_library(AMDiS::compositeFEM ALIAS amdis_compositeFEM)

target_compile_definitions(amdis_compositeFEM PUBLIC
    HAVE_COMPOSITE_FEM=1)
target_include_directories(amdis_compositeFEM PRIVATE ${COMPOSITE_SOURCE_DIR})
target_link_libraries(amdis_compositeFEM amdis)

# specify how to install this target:
# -----------------------------------

file(GLOB COMPOSITE_HEADERS "${COMPOSITE_SOURCE_DIR}/*.h*")
install(FILES ${COMPOSITE_HEADERS} DESTINATION include/amdis/compositeFEM)
install(TARGETS amdis_compositeFEM DESTINATION lib/amdis/ )