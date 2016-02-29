cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(AMDiS VERSION 2.0 LANGUAGES CXX)

if (POLICY CMP0060)
    cmake_policy(SET CMP0060 NEW)
endif ()

set(BASE_DIR ${AMDiS_SOURCE_DIR})
set(SOURCE_DIR ${BASE_DIR}/src)
list(APPEND CMAKE_MODULE_PATH ${BASE_DIR}/cmake)

if (NOT MSVC AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")
endif ()

# ------------------------------------------------------------------------------

#some options to control the build process
#option(ENABLE_PARALLEL_DOMAIN "Use parallel domain decomposition" false)
option(ENABLE_UMFPACK "Enable support for the UMFPACK solver" OFF)

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
    
target_compile_definitions(amdis_base INTERFACE AMDIS_FIXED_SIZE=1)
#target_compile_definitions(amdis_base INTERFACE DOW=2)

# ------------------------------------------------------------------------------

include(find_mtl4)  # -> target mtl4
install(FILES ${BASE_DIR}/cmake/find_mtl4.cmake DESTINATION share/amdis/)

include(find_umfpack)

include(find_boost) # -> target boost
install(FILES ${BASE_DIR}/cmake/find_boost.cmake DESTINATION share/amdis/)

# ------------------------------------------------------------------------------

add_library(amdis
    ${SOURCE_DIR}/AdaptBase.cpp
    ${SOURCE_DIR}/AdaptInfo.cpp
    ${SOURCE_DIR}/AdaptInstationary.cpp
    ${SOURCE_DIR}/AdaptStationary.cpp
    ${SOURCE_DIR}/AMDiS.cpp
    ${SOURCE_DIR}/AMDiS_base.cpp
    ${SOURCE_DIR}/Assembler.cpp
    ${SOURCE_DIR}/BasisFunction.cpp
    ${SOURCE_DIR}/Boundary.cpp
    ${SOURCE_DIR}/BoundaryManager.cpp
    ${SOURCE_DIR}/BoundaryObject.cpp
    ${SOURCE_DIR}/CoarseningManager.cpp
    ${SOURCE_DIR}/CoarseningManager1d.cpp
    ${SOURCE_DIR}/CoarseningManager2d.cpp
    ${SOURCE_DIR}/CoarseningManager3d.cpp
    ${SOURCE_DIR}/ComponentTraverseInfo.cpp
    ${SOURCE_DIR}/CouplingIterationInterface.cpp
    ${SOURCE_DIR}/CreatorMap.cpp
    ${SOURCE_DIR}/DOFAdmin.cpp
    #${SOURCE_DIR}/DOFIndexed.cpp
    ${SOURCE_DIR}/DOFMatrix.cpp
    ${SOURCE_DIR}/DOFVector.cpp
    ${SOURCE_DIR}/DOFVectorBase.cpp
    ${SOURCE_DIR}/DirichletBC.cpp
#     ${SOURCE_DIR}/DualTraverse.cpp
    ${SOURCE_DIR}/ElInfo.cpp
    ${SOURCE_DIR}/ElInfo1d.cpp
    ${SOURCE_DIR}/ElInfo2d.cpp
    ${SOURCE_DIR}/ElInfo3d.cpp
    ${SOURCE_DIR}/ElInfoStack.cpp
    ${SOURCE_DIR}/Element.cpp
    ${SOURCE_DIR}/ElementData.cpp
    ${SOURCE_DIR}/ElementDofIterator.cpp
    ${SOURCE_DIR}/ElementRegion_ED.cpp
    ${SOURCE_DIR}/FiniteElemSpace.cpp
    ${SOURCE_DIR}/FirstOrderAssembler.cpp
    ${SOURCE_DIR}/FirstOrderTerm.cpp
    ${SOURCE_DIR}/FixVec.cpp
    ${SOURCE_DIR}/Global.cpp
    ${SOURCE_DIR}/Initfile.cpp
    ${SOURCE_DIR}/Lagrange.cpp
    ${SOURCE_DIR}/LeafData.cpp
    ${SOURCE_DIR}/Line.cpp
    ${SOURCE_DIR}/Log.cpp
    ${SOURCE_DIR}/MacroElement.cpp
    ${SOURCE_DIR}/Marker.cpp
    ${SOURCE_DIR}/Mesh.cpp
    ${SOURCE_DIR}/MeshStructure.cpp
    ${SOURCE_DIR}/Operator.cpp
    ${SOURCE_DIR}/OperatorTerm.cpp
    ${SOURCE_DIR}/Parametric.cpp
    ${SOURCE_DIR}/PeriodicBC.cpp
    ${SOURCE_DIR}/ProblemInstat.cpp
    ${SOURCE_DIR}/ProblemStat.cpp
    ${SOURCE_DIR}/Projection.cpp
    ${SOURCE_DIR}/QPInfo.cpp
    ${SOURCE_DIR}/QPsiPhi.cpp
    ${SOURCE_DIR}/Quadrature.cpp
    ${SOURCE_DIR}/RCNeighbourList.cpp
    ${SOURCE_DIR}/RefinementManager.cpp
    ${SOURCE_DIR}/RefinementManager1d.cpp
    ${SOURCE_DIR}/RefinementManager2d.cpp
    ${SOURCE_DIR}/RefinementManager3d.cpp
    ${SOURCE_DIR}/RobinBC.cpp
    ${SOURCE_DIR}/ScalableQuadrature.cpp
    ${SOURCE_DIR}/SecondOrderAssembler.cpp
    ${SOURCE_DIR}/SecondOrderTerm.cpp
    ${SOURCE_DIR}/StandardProblemIteration.cpp
    ${SOURCE_DIR}/SubAssembler.cpp
    ${SOURCE_DIR}/SubElInfo.cpp
    ${SOURCE_DIR}/SubQuadrature.cpp
    ${SOURCE_DIR}/SurfaceQuadrature.cpp
    ${SOURCE_DIR}/SurfaceRegion_ED.cpp
    ${SOURCE_DIR}/SystemVector.cpp
    ${SOURCE_DIR}/Tetrahedron.cpp
    ${SOURCE_DIR}/Timer.cpp
    ${SOURCE_DIR}/Traverse.cpp
    ${SOURCE_DIR}/Triangle.cpp
    ${SOURCE_DIR}/VertexVector.cpp
    ${SOURCE_DIR}/ZeroOrderAssembler.cpp
    #${SOURCE_DIR}/ZeroOrderTerm.cpp
    ${SOURCE_DIR}/est/Estimator.cpp
    #${SOURCE_DIR}/est/RecoveryEstimator.cpp
    ${SOURCE_DIR}/est/ResidualEstimator.cpp
    ${SOURCE_DIR}/est/SimpleResidualEstimator.cpp
    #${SOURCE_DIR}/io/ArhReader.cpp
    #${SOURCE_DIR}/io/detail/ArhReader.cpp
    #${SOURCE_DIR}/io/Arh2Reader.cpp
    #${SOURCE_DIR}/io/Arh2Writer.cpp
    ${SOURCE_DIR}/io/Arh3Reader.cpp
    ${SOURCE_DIR}/io/Arh3Writer.cpp
    #${SOURCE_DIR}/io/detail/Arh2Reader.cpp
    #${SOURCE_DIR}/io/detail/ArhWriter.cpp
    #${SOURCE_DIR}/io/detail/Arh2Writer.cpp
    ${SOURCE_DIR}/io/detail/Arh3Reader.cpp
    ${SOURCE_DIR}/io/detail/Arh3Writer.cpp
    ${SOURCE_DIR}/io/DofWriter.cpp
    ${SOURCE_DIR}/io/ElementFileWriter.cpp
    ${SOURCE_DIR}/io/FileWriterInterface.cpp
    ${SOURCE_DIR}/io/FileWriter.cpp
    ${SOURCE_DIR}/io/GNUPlotWriter.cpp
    ${SOURCE_DIR}/io/MacroInfo.cpp
    ${SOURCE_DIR}/io/MacroReader.cpp
    ${SOURCE_DIR}/io/MacroWriter.cpp
    ${SOURCE_DIR}/io/PngReader.cpp
    ${SOURCE_DIR}/io/PngWriter.cpp
    ${SOURCE_DIR}/io/PovrayWriter.cpp
    ${SOURCE_DIR}/io/Spreadsheet.cpp
    ${SOURCE_DIR}/io/ValueReader.cpp
    ${SOURCE_DIR}/io/ValueWriter.cpp
    ${SOURCE_DIR}/io/VtkWriter.cpp
    ${SOURCE_DIR}/io/VtkVectorWriter.cpp
    ${SOURCE_DIR}/io/detail/VtkWriter.cpp
    #${SOURCE_DIR}/nonlin/ProblemNonLin.cpp
    ${SOURCE_DIR}/solver/LinearSolverInterface.cpp
    ${SOURCE_DIR}/solver/SolverMatrix.cpp
    ${SOURCE_DIR}/solver/Mapper.cpp
    #${SOURCE_DIR}/time/RosenbrockAdaptInstationary.cpp
    #${SOURCE_DIR}/time/RosenbrockMethod.cpp
    #${SOURCE_DIR}/time/RosenbrockStationary.cpp
)

add_library(amdis_debug INTERFACE)
target_sources(amdis_debug INTERFACE
    #${SOURCE_DIR}/Debug.cpp
    ${SOURCE_DIR}/ProblemStatDbg.cpp
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
