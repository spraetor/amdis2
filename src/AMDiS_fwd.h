/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file AMDiS_fwd.h */


#ifndef AMDIS_AMDIS_FWD_INCLUDE
#define AMDIS_AMDIS_FWD_INCLUDE

#include <boost/numeric/mtl/mtl.hpp>
#include "OpenMP.h"
#include "Config.h"

namespace AMDiS {

  class AdaptInfo;
  class AdaptStationary;
  class Assembler;
  class BasisFunction; 
  struct BlockMTLMatrix;
  class BoundaryManager;
//   class CGSolver;
  class CoarseningManager;
  class CouplingIterationInterface;
//   class DiagonalPreconditioner; 
  class DOFAdmin; 
  class DOFContainer;
  class DOFIndexedBase;
  class DOFMatrix;
  class DOFVectorDOF;
  class Element;
  class ElementDofIterator;
  class ElInfo; 
  class ElMatrixInfo;
  class Estimator;
  class FastQuadrature;
  class FillInfo; 
  class FileWriterInterface;
  class FiniteElemSpace; 
  class Flag;
//   class IdentityPreconditioner; 
//  class Parallel::InteriorBoundary;
  class LeafDataPeriodic;
  class LevelAdmin;
  class Line;
  class MacroElement;
  class MacroInfo;
  class Marker;
  class Mesh; 
  class MeshStructure;
  class LinearSolverInterface;
  class Operator;
  class OperatorTerm;
  class Parametric;
  class PeriodicBC;
  class ProblemInstat;
  class ProblemInstatScal;
  class ProblemInstatVec;
  class ProblemIterationInterface;
  class ProblemStatBase;
  class ProblemStatSeq;
  class ProblemTimeInterface;
  class Projection;
  class PreconditionerScal;
  class Quadrature; 
  class Q00PsiPhi;
  class Q0Psi;
  class Q10PsiPhi;
  class Q01PsiPhi;
  class Q1Psi;
  class RCNeighbourList;
  class RefinementManager;
  class RosenbrockAdaptInstationary;
  class RosenbrockStationary;
  class RobinBC;
  class SubElInfo;
  class SurfaceOperator;
  class SMIAdapter;
  class SystemVector;
  class Tetrahedron;
  class TraverseStack;
  class Triangle;
  class VertexInfo;
  class VertexVector;
  
  namespace io {
    class ElementFileWriter;
    class GNUPlotWriter;
    class MacroReader;
    struct MacroWriter;
    class PngWriter;
    class PovrayWriter;
    class Spreadsheet;
  }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
  namespace Parallel {
    class ElementObjectDatabase;
    class FeSpaceDofMap;
    class MatrixNnzStructure;
    class MeshLevelData;
    class PetscSolver;
    class PetscSolverFeti;
    class PetscSolverFetiDebug;
    class PetscSolverSchur;
    class PetscSolverGlobalMatrix;
    class PetscSolverGlobalBlockMatrix;
    class PetscSolverCahnHilliard;
    class PetscSolverNavierStokes;
    class ParallelDofMapping;
    class ParallelProblemStat;
    class ParallelSolver;
  }
#endif

  struct BoundaryObject;
  struct AtomicBoundary;
  
  template<typename T>                                 class DOFIndexed;
  template<typename T>                                 class DOFVectorBase;
  template<typename T>                                 class DOFVector;
  template<typename T>                                 class DimVec;
  template<typename T>                                 class DimMat;
  template<typename T>                                 class DirichletBC;
//   template<typename ITLSolver>                         class ITL_LinearSolverInterface;
  template<typename T, typename MatT, typename VecT >  class ITL_Preconditioner;
  template<typename T>                                 class Matrix;
  template<typename T>                                 class MatrixOfFixVecs;
  template<typename T>                                 class MatVecMultiplier;
  template<typename T>                                 class SolverMatrix;
  template<typename T>                                 class Vector;
  template<typename T>                                 class WorldVector;
  template<typename T>                                 class WorldMatrix;
  template<typename T>                                 class VectorOfFixVecs;
  
  namespace detail {
    template <class P> class CouplingProblemStat;
  }

  // some typedefs to introduce matrix and vector-types
  // ---------------------------------------------------------------------------
  
  template <class T>
  using DenseVector = mtl::dense_vector<T>;
  
  template <class T>
  using DenseMatrix = mtl::dense2D<T>;
  
  template <class T>
  using SparseMatrix = mtl::compressed2D<T>;

  // TODO: use templated ElementVector
  using ElementVector = DenseVector<double>;
  
  // TODO: use templated ElementMatrix;
  using ElementMatrix = DenseMatrix<double>;
  
  

} // namespace AMDiS

#endif // AMDIS_AMDIS_FWD_INCLUDE
