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


/** \file HypreSolver.h */


#ifndef AMDIS_HYPRE_SOLVER_H
#define AMDIS_HYPRE_SOLVER_H

#ifdef MTL_HAS_HYPRE

#include "solver/LinearSolver.h"
#include "solver/ITL_Preconditioner.h"
#include "MTL4Types.h"
#include "solver/itl/hypre.hpp"

#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/mtl/mtl.hpp>

namespace AMDiS {

  struct Hypre_Runner : public RunnerBase< MTLTypes::MTLMatrix, MTLTypes::MTLVector >
  {       
    typedef MTLTypes::MTLMatrix                   MatrixType;
    typedef MTLTypes::MTLVector                   VectorType;
    typedef RunnerBase< MatrixType, VectorType >  super;
    
    /** Interface to the HYPRE BoomerAMG solver [...]
     * Parameters provided by AMDiS:
     * 
     * [solver]->cycle mode: 
     * 	1...V-cycle
     *	2...W-cycle
     * 
     * [solver]->interpolation type:
     *  0...classical modified interpolation
     *	1...LS interpolation (for use with GSMG)
     * 	2...classical modified interpolation for hyperbolic PDEs
     * 	3...direct interpolation (with separation of weights)
     * 	4...multipass interpolation
     * 	5...multipass interpolation (with separation of weights)
     * 	6...extended+i interpolation
     * 	7...extended+i (if no common C neighbor) interpolation
     * 	8...standard interpolation
     * 	9...standard interpolation (with separation of weights)
     * 	10..classical block interpolation (for use with nodal systems version only)
     * 	11..classical block interpolation (for use with nodal systems version only)
     * 		with diagonalized diagonal blocks
     * 	12..FF interpolation
     * 	13..FF1 interpolation
     * 	14..extended interpolation
     * 
     *  [solver]->info: 
     * 	0...no printout (default) 
     * 	1...print setup information 
     * 	2...print solve information 
     * 	3...print both setup and solve information
     * 
     *  [solver]->relax type:
     * 	0...Jacobi
     * 	1...Gauss-Seidel, sequential (very slow!)
     * 	2...Gauss-Seidel, interior points in parallel, boundary sequential (slow!)
     * 	3...hybrid Gauss-Seidel or SOR, forward solve
     * 	4...hybrid Gauss-Seidel or SOR, backward solve
     * 	5...hybrid chaotic Gauss-Seidel (works only with OpenMP)
     * 	6...hybrid symmetric Gauss-Seidel or SSOR
     * 	9...Gaussian elimination (only on coarsest level)
     * */
    Hypre_Runner(LinearSolverInterface* oemPtr)
      : oem(*oemPtr),
	useTransposed(false),
	solverCreated(false)
    { 
      int cycleMode = -1, interpolation = -1, relaxation = -1;
      Parameters::get(oem.getName() + "->cycle mode", cycleMode);
      Parameters::get(oem.getName() + "->interpolation type", interpolation);
      Parameters::get(oem.getName() + "->relax type", relaxation);
            
      config.maxIter(oem.getMaxIterations());
      config.interpolationType(interpolation);
      config.relaxType(relaxation);
      config.cycle(cycleMode);
      config.tolerance(oem.getRelative());
      config.printLevel(oem.getInfo());
    }
    
    ~Hypre_Runner()
    {
      exit();
    }      
    
    /// Implementation of \ref RunnerBase::init()
    void init(const SolverMatrix<Matrix<DOFMatrix*> >& A, const MatrixType& mtlMatrix) override
    {
      setTransposed(typename MatrixType::orientation());
      // TODO: copy matrix directly from DOFMatrix to HYPRE matrix (?)
      hypreMatrix.init(mtlMatrix);
      HYPRE_IJMatrixGetObject(hypreMatrix, (void**) &matrix);
      HYPRE_BoomerAMGCreate(&solver);
      
      mtl::dense_vector<double> swap(1, 0.0);
      mtl::HypreParVector x(swap);
      HYPRE_BoomerAMGSetup(solver, matrix, x, x);
      
      solverCreated = true;
    }
        
    /// Implementation of \ref RunnerBase::solve()
    int solve(const MatrixType& A , VectorType& mtlX, const VectorType& mtlB) override
    {      
      mtl::HypreParVector x(mtlX);
      mtl::HypreParVector b(mtlB);
      config(solver);
      int error = 0;
      if(useTransposed)
	error = HYPRE_BoomerAMGSolveT(solver, matrix, b, x);
      else
	error = HYPRE_BoomerAMGSolve(solver, matrix, b, x);
      mtl::convert(x.getHypreVector(), mtlX);
      
      int num_iter = 0;
      HYPRE_BoomerAMGGetNumIterations(solver, &num_iter);
      oem.setIterations(num_iter);
      
      double rel_resid = 0.0;
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &rel_resid);
      oem.setRelativeResidual(rel_resid);
      
      oem.setErrorCode(error);
      return error;
    }
    
    /// Implementation of \ref RunnerInterface::exit()
    void exit()
    {
      if (solverCreated)
	HYPRE_BoomerAMGDestroy(solver);
      solverCreated = false;
    }
    
  private:
    inline void setTransposed(mtl::row_major)
    { 
      useTransposed = false;
    }
    
    inline void setTransposed(mtl::col_major)
    { 
      useTransposed = true;
    }
    
  protected:
    LinearSolverInterface& oem;
    
  private:
    HYPRE_Solver solver;
    HYPRE_ParCSRMatrix matrix;
    mtl::HypreMatrix hypreMatrix;
    
    mtl::AMGConfigurator config;
    
    bool useTransposed;
    bool solverCreated;
  };
  

  /**
   * \ingroup Solver
   * \class AMDiS::HypreSolver
   * \brief
   * Wrapper for the external HYPRE-AMG solver
   */
  typedef LinearSolver< MTLTypes::MTLMatrix, MTLTypes::MTLVector, Hypre_Runner > HypreSolver;

} // namespace AMDiS

#endif // MTL_HAS_HYPRE

#endif // AMDIS_HYPRE_SOLVER_H
