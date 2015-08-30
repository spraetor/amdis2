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


/** \file PITL_Solver.h */

#ifndef AMDIS_PITL_SOLVER_H
#define AMDIS_PITL_SOLVER_H

#ifdef HAVE_PARALLEL_MTL4

#include "solver/ITL_Solver.h"
#include "solver/LinearSolver.h"
#include "MTL4Types.h"

namespace AMDiS
{
  namespace Parallel
  {
    using namespace MTLTypes;

    template<typename SolverType>
    struct PITL_Solver : LinearSolver<PMTLMatrix, PMTLVector, ITL_Runner<SolverType, PMTLMatrix, PMTLVector>>
    {
      PITL_Solver(std::string name)
        : LinearSolver<PMTLMatrix, PMTLVector, ITL_Runner<SolverType, PMTLMatrix, PMTLVector>>(name) {}
    };

    typedef PITL_Solver<cg_solver_type> 	P_CGSolver;
    typedef PITL_Solver<cgs_solver_type> 	P_CGSSolver;
    //     typedef PITL_Solver< bicg_solver_type > 	P_BiCGSolver;
    typedef PITL_Solver<bicgstab_type> 	P_BiCGStabSolver;
    typedef PITL_Solver<bicgstab2_type> 	P_BiCGStab2Solver;
    typedef PITL_Solver<qmr_solver_type> 	P_QMRSolver;
    typedef PITL_Solver<tfqmr_solver_type> 	P_TFQMRSolver;
    typedef PITL_Solver<bicgstab_ell_type> 	P_BiCGStabEllSolver;
    typedef PITL_Solver<gmres_type> 		P_GMResSolver;
    typedef PITL_Solver<idr_s_type> 		P_IDRsSolver;
    typedef PITL_Solver<minres_solver_type> 	P_MinResSolver;
    typedef PITL_Solver<gcr_type> 		P_GcrSolver;
    typedef PITL_Solver<fgmres_type> 		P_FGMResSolver;
    typedef PITL_Solver<preonly_type> 	P_PreOnly;

    typedef ITL_PreconditionerBase<PMTLMatrix, PMTLVector> ParallelPreconditioner;
    typedef CreatorInterface<ParallelPreconditioner> ParallelPreconditionCreator;

    typedef ITL_Preconditioner<itl::pc::diagonal<PMTLMatrix>, PMTLMatrix, PMTLVector> 	P_DiagonalPreconditioner;
    typedef ITL_Preconditioner<itl::pc::identity<PMTLMatrix>, PMTLMatrix, PMTLVector> 	P_IdentityPreconditioner;
    typedef ITL_Preconditioner<itl::pc::ilu_0<PMTLMatrix>, PMTLMatrix, PMTLVector> 	P_ILUPreconditioner;
    typedef ITL_Preconditioner<itl::pc::ic_0<PMTLMatrix>, PMTLMatrix, PMTLVector> 	P_ICPreconditioner;

  } // end namespace Parallel

  template <>
  struct Collection<MTLTypes::PMTLMatrix>
  {
    typedef mtl::matrix::inserter<MTLTypes::PMTLMatrix , mtl::update_plus<MTLTypes::PMTLMatrix::value_type>> Inserter;
  };

  template <>
  struct Collection<MTLTypes::PMTLVector>
  {
    typedef mtl::vector::inserter<MTLTypes::PMTLVector, mtl::update_plus<MTLTypes::PMTLVector::value_type>> Inserter;
    typedef MTLTypes::PMTLMatrix PreconditionMatrix;
  };

} // namespace AMDiS

#endif // HAVE_PARALLEL_MTL4

#endif // AMDIS_PITL_SOLVER

