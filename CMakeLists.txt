cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(AMDiS VERSION 2.0 LANGUAGES CXX)

set(BASE_DIR ${AMDiS_SOURCE_DIR})
set(SOURCE_DIR ${BASE_DIR}/src)
list(APPEND CMAKE_MODULE_PATH ${BASE_DIR}/cmake)

if (NOT MSVC AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")
endif ()

# ------------------------------------------------------------------------------

#some options to control the build process
#option(ENABLE_PARALLEL_DOMAIN "Use parallel domain decomposition" false)
option(ENABLE_UMFPACK "Enable support for the UMFPACK solver" ON)

if (MSVC)
    option(BUILD_SHARED_LIBS "Build all libraries as shared or static, default: static" OFF)
    set(BUILD_SHARED_LIBS OFF)
    message(WARNING "Currently we can only build static libraries with Visual Studio")
else (MSVC)
    option(BUILD_SHARED_LIBS "Build all libraries as shared or static, default: shared" ON)
endif (MSVC)
mark_as_advanced(BUILD_SHARED_LIBS)


# ------------------------------------------------------------------------------

# define an interface-library for flags and compile-definitions
add_library(amdis_base INTERFACE)
add_library(AMDiS::base ALIAS amdis_base)

target_include_directories(amdis_base INTERFACE ${SOURCE_DIR})

target_compile_features(amdis_base INTERFACE
    cxx_alias_templates
    cxx_auto_type
    cxx_lambdas
    cxx_override
    cxx_range_for
    cxx_static_assert
    cxx_variadic_templates)

# ------------------------------------------------------------------------------

include(find_mtl4)  # -> target mtl4
install(FILES ${BASE_DIR}/cmake3/find_mtl4.cmake DESTINATION share/amdis/)

include(find_umfpack)

include(find_boost) # -> target boost
install(FILES ${BASE_DIR}/cmake3/find_boost.cmake DESTINATION share/amdis/)

# ------------------------------------------------------------------------------

add_library(amdis0
    ${SOURCE_DIR}/AdaptBase.cc)

add_library(amdis
    ${SOURCE_DIR}/AdaptBase.cc
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
    #${SOURCE_DIR}/DOFIndexed.cc
    ${SOURCE_DIR}/DOFMatrix.cc
    ${SOURCE_DIR}/DOFVector.cc
    ${SOURCE_DIR}/DOFVectorBase.cc
    ${SOURCE_DIR}/DirichletBC.cc
    #${SOURCE_DIR}/DualTraverse.cc
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
    #${SOURCE_DIR}/ZeroOrderTerm.cc
    ${SOURCE_DIR}/est/Estimator.cc
    #${SOURCE_DIR}/est/RecoveryEstimator.cc
    ${SOURCE_DIR}/est/ResidualEstimator.cc
    ${SOURCE_DIR}/est/SimpleResidualEstimator.cc
    #${SOURCE_DIR}/io/ArhReader.cc
    #${SOURCE_DIR}/io/detail/ArhReader.cc
    #${SOURCE_DIR}/io/Arh2Reader.cc
    #${SOURCE_DIR}/io/Arh2Writer.cc
    ${SOURCE_DIR}/io/Arh3Reader.cc
    ${SOURCE_DIR}/io/Arh3Writer.cc
    #${SOURCE_DIR}/io/detail/Arh2Reader.cc
    #${SOURCE_DIR}/io/detail/ArhWriter.cc
    #${SOURCE_DIR}/io/detail/Arh2Writer.cc
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
    #${SOURCE_DIR}/nonlin/ProblemNonLin.cc
    ${SOURCE_DIR}/solver/LinearSolverInterface.cc
    ${SOURCE_DIR}/solver/SolverMatrix.cc
    ${SOURCE_DIR}/solver/Mapper.cc
    #${SOURCE_DIR}/time/RosenbrockAdaptInstationary.cc
    #${SOURCE_DIR}/time/RosenbrockMethod.cc
    #${SOURCE_DIR}/time/RosenbrockStationary.cc
)

add_library(amdis_debug INTERFACE)
target_sources(amdis_debug INTERFACE
    #${SOURCE_DIR}/Debug.cc
    ${SOURCE_DIR}/ProblemStatDbg.cc
)

#include(amdis_parallel)     # adds sources to amdis
#include(amdis_extensions)   # -> target AMDiS::extensions
include(muparser)           # -> target muparser

# ------------------------------------------------------------------------------

if (ENABLE_COMPRESSION)
    target_compile_definitions(amdis PRIVATE HAVE_COMPRESSION=1)
endif (ENABLE_COMPRESSION)

if (ENABLE_UMFPACK)
    target_compile_definitions(amdis PRIVATE HAVE_UMFPACK=1)
endif (ENABLE_UMFPACK)

# if (ENABLE_REINIT)
#     include(amdis_reinit)       # -> target AMDiS::reinit
# endif (ENABLE_REINIT)
#
# if (ENABLE_COMPOSITE_FEM)
#     include(amdis_compositeFEM) # -> target AMDiS::compositeFEM
# endif (ENABLE_COMPOSITE_FEM)

# ------------------------------------------------------------------------------

if (MSVC)
    target_compile_definitions(amdis PRIVATE
        _SCL_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_WARNINGS)
endif (MSVC)

target_link_libraries(amdis0 AMDiS::base)

target_link_libraries(amdis AMDiS::base muparser
    $<$<NOT:$<CONFIG:Release>>:amdis_debug>)

# specify how to install this target:
# -----------------------------------

set(INSTALL_SUBDIRS . config nonlin est expressions
    matrix_vector matrix_vector/expr operations traits
	utility solver solver/details solver/itl io io/detail)

foreach (SUBDIR ${INSTALL_SUBDIRS})
    file(GLOB HEADERS "${SOURCE_DIR}/${SUBDIR}/*.h*")
    install(FILES ${HEADERS} DESTINATION include/amdis/${SUBDIR}/)
endforeach ()

install(TARGETS amdis DESTINATION lib/amdis/ )

# generate configuration file:
# ----------------------------

configure_file(${BASE_DIR}/cmake/AMDISConfig.cmake.in
    ${AMDiS_BINARY_DIR}/AMDISConfig.cmake
    @ONLY
)

install(FILES ${AMDiS_BINARY_DIR}/AMDISConfig.cmake DESTINATION share/amdis/)
install(FILES ${BASE_DIR}/cmake/AMDISUse.cmake DESTINATION share/amdis/)
install(FILES ${BASE_DIR}/cmake/test_compile_features.cmake DESTINATION share/amdis/)
install(FILES ${BASE_DIR}/cmake/dummy.cc DESTINATION share/amdis/)
