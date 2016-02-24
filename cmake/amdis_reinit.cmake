# specify the target and requirements for the reinit library

SET(COMPOSITE_SOURCE_DIR ${SOURCE_DIR}/reinit)
set(REINIT_SOURCE_DIR ${SOURCE_DIR}/reinit)
file(GLOB REINIT_SRC ${REINIT_SOURCE_DIR}/*.cc)
add_library(amdis_reinit ${REINIT_SRC})
add_library(AMDiS::reinit ALIAS amdis_reinit)

target_compile_definitions(amdis_reinit PUBLIC
    HAVE_REINIT=1)
target_include_directories(amdis_reinit PUBLIC ${REINIT_SOURCE_DIR})
target_link_libraries(amdis_reinit amdis)

# specify how to install this target:
# -----------------------------------

file(GLOB RINIT_HEADERS "${REINIT_SOURCE_DIR}/*.h*")
install(FILES ${RINIT_HEADERS} DESTINATION include/amdis/reinit)
install(TARGETS amdis_reinit DESTINATION lib/amdis/ )