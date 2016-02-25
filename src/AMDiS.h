/** \file AMDiS.h */

#pragma once

// std c++ headers
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <cstdint>

// AMDiS includes
#include "AMDiS.hpp"
#include "MTL4Types.hpp"
#include "AdaptInfo.hpp"
#include "AdaptInstationary.hpp"
#include "AdaptStationary.hpp"
#include "Assembler.hpp"
#include "BasisFunction.hpp"
#include "Boundary.hpp"
#include "Boundary.hpp"
#include "BoundaryCondition.hpp"
#include "BoundaryManager.hpp"
#include "CoarseningManager.hpp"
#include "CoarseningManager1d.hpp"
#include "CoarseningManager2d.hpp"
#include "CoarseningManager3d.hpp"
#include "CouplingTimeInterface.hpp"
#include "CouplingIterationInterface.hpp"
#include "CreatorInterface.hpp"
#include "CreatorMap.hpp"
#include "Debug.hpp"
#include "DOFAdmin.hpp"
#include "DOFContainer.hpp"
#include "DOFIndexed.hpp"
#include "DOFIterator.hpp"
#include "DOFMatrix.hpp"
#include "DOFVector.hpp"
#include "DOFVectorOperations.hpp"
#include "DirichletBC.hpp"
// #include "DualTraverse.hpp"
#include "ElInfo.hpp"
#include "ElInfo1d.hpp"
#include "ElInfo2d.hpp"
#include "ElInfo3d.hpp"
#include "Element.hpp"
#include "ElementDofIterator.hpp"
#include "FiniteElemSpace.hpp"
#include "FirstOrderTerm.hpp"
#include "FixVec.hpp"
#include "Flag.hpp"
#include "Global.hpp"
#include "Initfile.hpp"
#include "Lagrange.hpp"
#include "LeafData.hpp"
#include "Line.hpp"
#include "MacroElement.hpp"
#include "Marker.hpp"
#include "MatrixVector.hpp"
#include "MatrixVectorOperations.hpp"
#include "Mesh.hpp"
#include "MeshStructure.hpp"
#include "ComponentTraverseInfo.hpp"
#include "Operator.hpp"
#include "OperatorTerm.hpp"
#include "Parametric.hpp"
#include "PeriodicBC.hpp"
#include "ProblemStat.hpp"
#include "ProblemInstat.hpp"
#include "ProblemTimeInterface.hpp"
#include "ProblemStatBase.hpp"
#include "SecondOrderTerm.hpp"
#include "StandardProblemIteration.hpp"
#include "Projection.hpp"
#include "QPsiPhi.hpp"
#include "Quadrature.hpp"
#include "RCNeighbourList.hpp"
#include "RefinementManager.hpp"
#include "RefinementManager1d.hpp"
#include "RefinementManager2d.hpp"
#include "RefinementManager3d.hpp"
#include "RobinBC.hpp"
#include "SurfaceOperator.hpp"
#include "SurfaceQuadrature.hpp"
#include "SystemVector.hpp"
#include "Tetrahedron.hpp"
#include "Traverse.hpp"
#include "Traits.hpp"
#include "Triangle.hpp"
#include "VertexVector.hpp"
#include "ZeroOrderTerm.hpp"

#include "est/Estimator.hpp"

// #include "io/ArhReader.hpp"
// #include "io/Arh2Reader.hpp"
#include "io/Arh3Reader.hpp"
// #include "io/ArhWriter.hpp"
// #include "io/Arh2Writer.hpp"
#include "io/Arh3Writer.hpp"
#include "io/DataCollector.hpp"
#include "io/FileWriter.hpp"
#include "io/GNUPlotWriter.hpp"
#include "io/GridWriter.hpp"
#include "io/MacroWriter.hpp"
#include "io/PngWriter.hpp"
#include "io/PovrayWriter.hpp"
#include "io/Spreadsheet.hpp"
#include "io/ValueReader.hpp"
#include "io/ValueWriter.hpp"
#include "io/VtkWriter.hpp"
#include "io/VtkVectorWriter.hpp"
//#include "io/VtkReader.h"
#include "io/Reader.hpp"
#include "io/Writer.hpp"

// #include "nonlin/ProblemNonLin.h"
// #include "nonlin/NonLinSolver.h"

#include "solver/ITL_Preconditioner.hpp"
#include "solver/ITL_Solver.hpp"
#include "solver/LinearSolverInterface.hpp"

// #include "time/RosenbrockAdaptInstationary.h"
// #include "time/RosenbrockStationary.h"


#if HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/InteriorBoundary.hpp"
#include "parallel/MpiHelper.hpp"
#include "parallel/ParallelDebug.hpp"
#include "parallel/StdMpi.hpp"
#include "parallel/ParallelProblemStat.hpp"

#if !HAVE_PARALLEL_MTL4
#include "parallel/PetscSolver.hpp"
#include "parallel/PetscSolverNavierStokes.hpp"
#endif

#endif

#if HAVE_PETSC
#include <petsc.h>
#endif
