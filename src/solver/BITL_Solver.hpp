/** \file BITL_Solver.h */

#pragma once

#include "MTL4Types.hpp"
#include "solver/BlockMTLMatrix.hpp"
#include "solver/LinearSolver.hpp"
#include "solver/ITL_Runner.hpp"
#include "solver/ITL_Solver.hpp"
#include "solver/itl/block_diagonal.hpp"

namespace AMDiS
{
  /**
   * \ingroup Solver
   *
   * \brief
   * Wrapper for MTL4 itl-solvers applied to block matrices.
   *
   */
  
  template <class SolverType>
  using BITL_Solver = 
    LinearSolver<BlockMTLMatrix, 
		 MTLTypes::MTLVector, 
		 ITL_Runner<SolverType, BlockMTLMatrix, MTLTypes::MTLVector>>;


  // using typedefs to specify the krylov-subspace algorithm for the solvers
  // ---------------------------------------------------------------------------

  using B_CGSolver       = BITL_Solver<cg_solver_type>;
  using B_CGSSolver      = BITL_Solver<cgs_solver_type>;
  using B_BiCGStabSolver = BITL_Solver<bicgstab_type>;
  using B_TFQMRSolver    = BITL_Solver<tfqmr_solver_type>;
  using B_GMResSolver    = BITL_Solver<gmres_type>;
  using B_MinResSolver   = BITL_Solver<minres_solver_type>;
  using B_GcrSolver      = BITL_Solver<gcr_type>;
  using B_FGMResSolver   = BITL_Solver<fgmres_type>;
  using B_PreOnly        = BITL_Solver<preonly_type>;

  //   using B_BiCGSolver = BITL_Solver< bicg_solver_type >;   	// uses adjoint(A)
  //   using B_BiCGStab2Solver = BITL_Solver< bicgstab2_type >;  	// ERROR
  //   using B_QMRSolver = BITL_Solver< qmr_solver_type >;    	// uses trans(A)
  //   using B_BiCGStabEllSolver = BITL_Solver< bicgstab_ell_type >; 	// ERROR
  //   using B_IDRsSolver = BITL_Solver< idr_s_type >;  		// ERROR


  // using typedefs to specify some available preconditioners for block systems
  // ---------------------------------------------------------------------------
  
  template <class Precon>
  using BITL_Preconditioner =
    ITL_Preconditioner<Precon, BlockMTLMatrix, MTLTypes::MTLVector>;

  using BlockDiagonalPreconditioner = 
    BITL_Preconditioner<itl::pc::diagonal<BlockMTLMatrix>>;
  using BlockIdentityPreconditioner = 
    BITL_Preconditioner<itl::pc::identity<BlockMTLMatrix>>;


  /// Function to initialize the creatormap for arbitrary matrix-types
  template <class MatrixType>
  void initCreatorMap(std::string backend)
  {
    using namespace MTLTypes;
    template <class SolverType>
    using Solver = 
      LinearSolver<MatrixType, MTLVector, ITL_Runner<SolverType, MatrixType, MTLVector>>;
    
    using MatCgSolver       = typename Solver<cg_solver_type>::Creator;
    using MatCgsSolver      = typename Solver<cgs_solver_type>::Creator;
    using MatBicgstabSolver = typename Solver<bicgstab_type>::Creator;
    using MatTfqmrSolver    = typename Solver<tfqmr_solver_type>::Creator;
    using MatGmresSolver    = typename Solver<gmres_type>::Creator;
    using MatMinresSolver   = typename Solver<minres_solver_type>::Creator;
    using MatGcrSolver      = typename Solver<gcr_type>::Creator;
    using MatPreonlySolver  = typename Solver<preonly_type>::Creator;
    using MatFgmresSolver   = typename Solver<fgmres_type>::Creator;

    CreatorMap<LinearSolverInterface>::addCreator(backend + "_cg",       new MatCgSolver);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_cgs",      new MatCgsSolver);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_bicgstab", new MatBicgstabSolver);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_tfqmr",    new MatTfqmrSolver);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_gmres",    new MatGmresSolver);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_minres",   new MatMinresSolver);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_gcr",      new MatGcrSolver);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_preconly", new MatPreonlySolver);
    CreatorMap<LinearSolverInterface>::addCreator(backend + "_fgmres",   new MatFgmresSolver);
  }

} // end namespace AMDiS
