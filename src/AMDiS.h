/** \file AMDiS.h */

#pragma once

#include "stdint.h"
#include "MTL4Types.h"
#include "AdaptInfo.h"
#include "AdaptInstationary.h"
#include "AdaptStationary.h"
#include "Assembler.h"
#include "BasisFunction.h"
#include "Boundary.h"
#include "Boundary.h"
#include "BoundaryCondition.h"
#include "BoundaryManager.h"
#include "CoarseningManager.h"
#include "CoarseningManager1d.h"
#include "CoarseningManager2d.h"
#include "CoarseningManager3d.h"
#include "CouplingTimeInterface.h"
#include "CouplingIterationInterface.h"
#include "CreatorInterface.h"
#include "CreatorMap.h"
#include "Debug.h"
#include "DOFAdmin.h"
#include "DOFContainer.h"
#include "DOFIndexed.h"
#include "DOFIterator.h"
#include "DOFMatrix.h"
#include "DOFVector.h"
#include "DirichletBC.h"
#include "DualTraverse.h"
#include "ElInfo.h"
#include "ElInfo1d.h"
#include "ElInfo2d.h"
#include "ElInfo3d.h"
#include "Element.h"
#include "ElementDofIterator.h"
#include "Error.h"
#include "FiniteElemSpace.h"
#include "FirstOrderTerm.h"
#include "FixVec.h"
#include "Flag.h"
#include "Global.h"
#include "Initfile.h"
#include "Lagrange.h"
#include "LeafData.h"
#include "Line.h"
#include "MacroElement.h"
#include "Marker.h"
#include "MatrixVector.h"
#include "MatrixVectorOperations.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "ComponentTraverseInfo.h"
#include "Operator.h"
#include "OperatorTerm.h"
#include "Parametric.h"
#include "PeriodicBC.h"
#include "ProblemStat.h"
#include "ProblemInstat.h"
#include "ProblemTimeInterface.h"
#include "ProblemInterpol.h"
#include "ProblemStatBase.h"
#include "SecondOrderTerm.h"
#include "StandardProblemIteration.h"
#include "Projection.h"
#include "QPsiPhi.h"
#include "Quadrature.h"
#include "RCNeighbourList.h"
#include "RefinementManager.h"
#include "RefinementManager1d.h"
#include "RefinementManager2d.h"
#include "RefinementManager3d.h"
#include "RobinBC.h"
#include "SurfaceOperator.h"
#include "SurfaceQuadrature.h"
#include "SystemVector.h"
#include "Tetrahedron.h"
#include "TimedObject.h"
#include "Traverse.h"
#include "Traits.h"
#include "Triangle.h"
#include "VertexVector.h"
#include "ZeroOrderTerm.h"

#include "est/Estimator.h"

#include "io/ArhReader.h"
#include "io/Arh2Reader.h"
#include "io/Arh3Reader.h"
#include "io/ArhWriter.h"
#include "io/Arh2Writer.h"
#include "io/Arh3Writer.h"
#include "io/DataCollector.h"
#include "io/FileWriter.h"
#include "io/GNUPlotWriter.h"
#include "io/GridWriter.h"
#include "io/MacroWriter.h"
#include "io/PngWriter.h"
#include "io/PovrayWriter.h"
#include "io/Spreadsheet.h"
#include "io/ValueReader.h"
#include "io/ValueWriter.h"
#include "io/VtkWriter.h"
#include "io/VtkVectorWriter.h"
#ifdef HAVE_EXTENSIONS
#include "io/VtkReader.h"
#endif
#include "io/Reader.h"
#include "io/Writer.h"

#include "nonlin/ProblemNonLin.h"
#include "nonlin/NonLinSolver.h"

#include "solver/ITL_Preconditioner.h"
#include "solver/ITL_Solver.h"
#include "solver/LinearSolverInterface.h"

#include "time/RosenbrockAdaptInstationary.h"
#include "time/RosenbrockStationary.h"


#if HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/InteriorBoundary.h"
#include "parallel/MpiHelper.h"
#include "parallel/ParallelDebug.h"
#include "parallel/StdMpi.h"
#include "parallel/ParallelProblemStat.h"

#if !HAVE_PARALLEL_MTL4
#include "parallel/PetscSolver.h"
#include "parallel/PetscSolverNavierStokes.h"
#endif

#endif

#if HAVE_PETSC
#include <petsc.h>
#endif

namespace AMDiS {

  void init(int argc, char **argv, std::string initFileName = "");

  void init(std::string initFileName);

  void finalize();
  
} // end namespace AMDiS
