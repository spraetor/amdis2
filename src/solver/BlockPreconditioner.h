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
 

/** \file BlockPreconditioner.h */

#ifndef AMDIS_BLOCK_PRECONDITIONER_H
#define AMDIS_BLOCK_PRECONDITIONER_H

#include "solver/ITL_Preconditioner.h"
#include "solver/SolverMatrix.h"

namespace AMDiS {

  /// Basis preconditioner structure for block-preconditioners
  template <class MatrixType, class VectorType = MTLTypes::MTLVector>
  struct BlockPreconditioner : public ITL_PreconditionerBase<MatrixType, VectorType>
  {    
    typedef BlockPreconditioner                             self;
    typedef ITL_PreconditionerBase<MatrixType, VectorType>  super;
    typedef super                                           precon_base;
        
    BlockPreconditioner() 
      : A(NULL), fullMatrix(NULL)
    { }
    
    virtual ~BlockPreconditioner() {}
    
    /// extract iranges from BlockMatrix to be used to extract sub-vectors and sub-matrices.
    void init(const SolverMatrix<Matrix<DOFMatrix*> >& A_, const MatrixType& fullMatrix_) override
    {
      A = &A_;
      fullMatrix = &fullMatrix_;
      
      BlockMapper mapper(A_);
      rows.resize(mapper.getNumComponents());
      int start = 0;
      for (int i = 0; i < mapper.getNumComponents(); i++) {
	mapper.setRow(i+1);
	int finish = mapper.row(0);
	rows[i].set(start, finish);
	start = finish;
      }
    }
  
    virtual void solve(const VectorType& b, VectorType& x) const
    { FUNCNAME("BlockPreconditioner::solve()");
      TEST_EXIT(false)("Must be implemented in derived classes!\n");
    }
    
    virtual void adjoint_solve(const VectorType& x, VectorType& y) const
    { FUNCNAME("BlockPreconditioner::adjoint_solve()");
      TEST_EXIT(false)("Must be implemented in derived classes!\n");
    }
    
    inline const mtl::irange& getRowRange(size_t i) const
    {
      return rows[i];
    }
    
    inline const mtl::irange& getColRange(size_t i) const
    {
      return rows[i];
    } 
    
    const DOFMatrix::base_matrix_type& getSubMatrix(size_t i, size_t j) const 
    {
      return A->getSubMatrix(i, j);
    }
    
    template <class SolverType, class RunnerType>
    static void createSubSolver(std::string param, SolverType*& solver, RunnerType*& runner, 
				std::string solverType = "0", std::string preconType = "no", 
				int max_iter = 100, double tol = 1.e-8)
    {
      // definition of standard-backends
#if defined HAVE_PARALLEL_PETSC
      std::string backend("p_petsc");
#elif defined HAVE_PARALLEL_MTL
      std::string backend("p_mtl");
#elif defined HAVE_PETSC
      std::string backend("petsc");
#else
      std::string backend("mtl");
#endif
    
      // === read backend-name ===
      std::string initFileStr = param + "->solver";    

      // === read solver-name ===
      std::string tmp = "0";
      Parameters::get(initFileStr, tmp);
      if (tmp == "0") {
	TEST_EXIT(solverType != "0")("You have to specify the parameter %s in the initfile!\n", initFileStr.c_str());
	
	Parameters::set(initFileStr, solverType);
	Parameters::set(initFileStr + "->backend", backend);
	Parameters::set(initFileStr + "->left precon", preconType);
	Parameters::set(initFileStr + "->max iteration", max_iter);
	Parameters::set(initFileStr + "->tolerance", tol);
      } else {
	solverType = tmp;
      }
      
      Parameters::get(initFileStr + "->backend", backend);
      
      if (backend != "0" && backend != "no" && backend != "")
	solverType = backend + "_" + solverType;
    
      LinearSolverCreator *solverCreator = 
	dynamic_cast<LinearSolverCreator*>(CreatorMap<LinearSolverInterface>::getCreator(solverType, initFileStr));
      TEST_EXIT(solverCreator)
	("No valid solver type found in parameter \"%s\"\n", initFileStr.c_str());
      solverCreator->setName(initFileStr);
      solver = dynamic_cast<SolverType*>(solverCreator->create());
      assert(solver != NULL);
      
      runner = dynamic_cast<RunnerType*>(solver->getRunner());
      assert(runner != NULL);
    }
    
    
  protected:
    const SolverMatrix<Matrix<DOFMatrix*> >* A;
    const MatrixType* fullMatrix;
    
    std::vector<mtl::irange> rows;
  };
  
} // end namespace AMDiS

#endif // AMDIS_BLOCK_PRECONDITIONER_H
