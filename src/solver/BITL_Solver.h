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


/** \file BITL_Solver.h */

#ifndef AMDIS_BITL_SOLVER_H
#define AMDIS_BITL_SOLVER_H

#include "solver/BlockMTLMatrix.h"

#include "solver/LinearSolver.h"
#include "solver/ITL_Runner.h"
#include "solver/ITL_Solver.h"
#include "MTL4Types.h"

#include "solver/itl/block_diagonal.hpp"

namespace AMDiS {
  

  /** 
   * \ingroup Solver
   * 
   * \brief 
   * Wrapper for MTL4 itl-solvers applied to block matrices. 
   * 
   */  
  template< typename SolverType >
  struct BITL_Solver : LinearSolver< BlockMTLMatrix, MTLTypes::MTLVector, ITL_Runner< SolverType, BlockMTLMatrix, MTLTypes::MTLVector > >
  {
    BITL_Solver(std::string name)
      : LinearSolver< BlockMTLMatrix, MTLTypes::MTLVector, ITL_Runner< SolverType, BlockMTLMatrix, MTLTypes::MTLVector > >(name) 
    { }
  };
  

  // ===========================================================================
  // using typedefs to specify the krylov-subspace algorithm for the solvers

  typedef BITL_Solver< cg_solver_type >    B_CGSolver;
  typedef BITL_Solver< cgs_solver_type >   B_CGSSolver;
  typedef BITL_Solver< bicgstab_type >     B_BiCGStabSolver;
  typedef BITL_Solver< tfqmr_solver_type > B_TFQMRSolver;
  typedef BITL_Solver< gmres_type >        B_GMResSolver;
  typedef BITL_Solver< minres_solver_type> B_MinResSolver;
  typedef BITL_Solver< gcr_type >          B_GcrSolver;
  typedef BITL_Solver< fgmres_type >       B_FGMResSolver;
  typedef BITL_Solver< preonly_type >      B_PreOnly;
  
  //   typedef BITL_Solver< bicg_solver_type >  B_BiCGSolver;   	// uses adjoint(A)
  //   typedef BITL_Solver< bicgstab2_type >    B_BiCGStab2Solver;  	// ERROR
  //   typedef BITL_Solver< qmr_solver_type >   B_QMRSolver;    	// uses trans(A)
  //   typedef BITL_Solver< bicgstab_ell_type > B_BiCGStabEllSolver; 	// ERROR
  //   typedef BITL_Solver< idr_s_type >        B_IDRsSolver;  		// ERROR
  
  
  // ===========================================================================
  // using typedefs to specify some available preconditioners for block systems

  typedef ITL_Preconditioner<itl::pc::diagonal<BlockMTLMatrix>, BlockMTLMatrix, MTLTypes::MTLVector > BlockDiagonalPreconditioner;
  typedef ITL_Preconditioner<itl::pc::identity<BlockMTLMatrix>, BlockMTLMatrix, MTLTypes::MTLVector > BlockIdentityPreconditioner;
  
  
  /// Function to initialize the creatormap for arbitrary matrix-types
  template<typename MatrixType>
  void initCreatorMap(std::string backend)
  {
    using namespace MTLTypes;
    typedef LinearSolver<MatrixType, MTLVector, ITL_Runner<cg_solver_type, MatrixType, MTLVector> >    MatCgSolver;
    typedef LinearSolver<MatrixType, MTLVector, ITL_Runner<cgs_solver_type, MatrixType, MTLVector> >   MatCgsSolver;
    typedef LinearSolver<MatrixType, MTLVector, ITL_Runner<bicgstab_type, MatrixType, MTLVector> >     MatBicgstabSolver;
    typedef LinearSolver<MatrixType, MTLVector, ITL_Runner<tfqmr_solver_type, MatrixType, MTLVector> > MatTfqmrSolver;
    typedef LinearSolver<MatrixType, MTLVector, ITL_Runner<gmres_type, MatrixType, MTLVector> >        MatGmresSolver;
    typedef LinearSolver<MatrixType, MTLVector, ITL_Runner<minres_solver_type, MatrixType, MTLVector> > MatMinresSolver;
    typedef LinearSolver<MatrixType, MTLVector, ITL_Runner<gcr_type, MatrixType, MTLVector> >          MatGcrSolver;
    typedef LinearSolver<MatrixType, MTLVector, ITL_Runner<preonly_type, MatrixType, MTLVector> >      MatPreonlySolver;
    typedef LinearSolver<MatrixType, MTLVector, ITL_Runner<fgmres_type, MatrixType, MTLVector> >       MatFgmresSolver;
    
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_cg", new typename MatCgSolver::Creator);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_cgs", new typename MatCgsSolver::Creator);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_bicgstab", new typename MatBicgstabSolver::Creator);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_tfqmr", new typename MatTfqmrSolver::Creator);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_gmres", new typename MatGmresSolver::Creator);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_minres", new typename MatMinresSolver::Creator);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_gcr", new typename MatGcrSolver::Creator);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_preconly", new typename MatPreonlySolver::Creator);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_fgmres", new typename MatFgmresSolver::Creator);
  }
  
} // end namespace AMDiS

#endif // AMDIS_BITL_SOLVER
 
