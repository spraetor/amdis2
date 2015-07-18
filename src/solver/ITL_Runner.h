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


/** \file ITL_Runner.h */


#ifndef AMDIS_ITL_RUNNER_H
#define AMDIS_ITL_RUNNER_H

#include "solver/LinearSolver.h"
#include "solver/ITL_Preconditioner.h"

#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/mtl/mtl.hpp>

namespace AMDiS {

  template< typename MatrixType, typename VectorType >
  struct PreconPair
  { 
    /// Pointer to the left preconditioner
    ITL_PreconditionerBase<MatrixType, VectorType>* l;

    /// Pointer to the right preconditioner
    ITL_PreconditionerBase<MatrixType, VectorType>* r;
    
    PreconPair() 
      : l(NULL), r(NULL) 
    { }
  };


  /** \ingroup Solver
   * 
   * \brief
   * Wrapper class for different MTL4 itl-solvers. These solvers
   * are parametrized by Matrix- and VectorType. 
   **/
  template< typename ITLSolver, typename MatrixType, typename VectorType >
  struct ITL_Runner : public RunnerBase< MatrixType, VectorType >
  {       
    typedef RunnerBase< MatrixType, VectorType > super;
    
    ITL_Runner(LinearSolverInterface* oemPtr)
      : oem(*oemPtr),
	solver(oem.getName())
    {
      setPrecon(preconPair);
    }
    
    ~ITL_Runner()
    {
      if (preconPair.l != NULL) {
	preconPair.l->exit();
	delete preconPair.l;
	preconPair.l = NULL;
      }

      if (preconPair.r != NULL) {
	preconPair.r->exit();
	delete preconPair.r;
	preconPair.r = NULL;
      }
    }
    
    
    /// Implementation of \ref RunnerBase::init()
    virtual void init(const SolverMatrix<Matrix<DOFMatrix*> >& A, const MatrixType& fullMatrix) override
    {
      preconPair.l->init(A, fullMatrix);
      preconPair.r->init(A, fullMatrix);
    }
    
    
    /// Implementation of \ref RunnerBase::solve()
    virtual int solve(const MatrixType& A , VectorType& x, const VectorType& b) override
    { FUNCNAME("ITL_Runner::solve()");
    
      TEST_EXIT(preconPair.l != NULL)("there is no left preconditioner\n");
      TEST_EXIT(preconPair.r != NULL)("there is no right preconditioner\n");
      
      typedef typename mtl::Collection<MatrixType>::value_type value_type;
            
      VectorType r(num_rows(b)); r = b;
      if (two_norm(x) != 0) {
	r = A * x ; 
	r -= b;
      }
      int error = 0;
      if (oem.getInfo() == 0) {
        // iteration that does not print residual per iteration
	itl::basic_iteration<value_type> 
	  iter(r, oem.getMaxIterations(), oem.getRelative(), oem.getTolerance());
	
	error = solver(A, x, b, *(preconPair.l), *(preconPair.r), iter);
	oem.setErrorCode(error);
	oem.setIterations(iter.iterations());
	oem.setResidual(iter.resid());
	oem.setRelativeResidual(iter.relresid());
      } else {
        // print information about the solution process
	itl::cyclic_iteration<value_type> 
	  iter(r, oem.getMaxIterations(), oem.getRelative(), oem.getTolerance(), 
	       oem.getPrint_cycle());
	
	error = solver(A, x, b, *(preconPair.l), *(preconPair.r), iter);
	oem.setErrorCode(error);
	oem.setIterations(iter.iterations());
	oem.setResidual(iter.resid());
	oem.setRelativeResidual(iter.relresid());
      }
      
      return error;
    }

    
    /// Implementation of \ref RunnerBase::adjoint_solve()
    virtual int adjoint_solve(const MatrixType& A , VectorType& x, const VectorType& b) override
    { FUNCNAME("ITL_Runner::adjoint_solve()");
    
      TEST_EXIT(preconPair.l != NULL)("there is no left preconditioner\n");
      TEST_EXIT(preconPair.r != NULL)("there is no right preconditioner\n");
            
#if 0
      typedef typename mtl::Collection<MatrixType>::value_type value_type;
      mtl::matrix::transposed_view<const MatrixType> B(A);
      VectorType r(B * x - b); 
      int error = 0;
      if (oem.getInfo() == 0) {
	itl::basic_iteration<value_type> 
	  iter(r, oem.getMaxIterations(), oem.getRelative(), oem.getTolerance());
	
	error = solver(B, x, b, *(preconPair.l), *(preconPair.r), iter);
	oem.setErrorCode(error);
	oem.setIterations(iter.iterations());
	oem.setResidual(iter.resid());
      } else {
	itl::cyclic_iteration<value_type> 
	  iter(r, oem.getMaxIterations(), oem.getRelative(), oem.getTolerance(), 
	       oem.getPrint_cycle());
	
	error = solver(B, x, b, *(preconPair.l), *(preconPair.r), iter);
	oem.setErrorCode(error);
	oem.setIterations(iter.iterations());
	oem.setResidual(iter.resid());
      }
      
      return error;
#endif
      return 1;
    }
      
      
    /// Implementation of \ref RunnerInterface::exit()
    virtual void exit() override
    {
      preconPair.l->exit();
      preconPair.r->exit();
    }
    
    
    /// Implementation of \ref RunnerInterface::getLeftPrecon()
    virtual PreconditionerInterface* getLeftPrecon() override
    {
      return preconPair.l;
    }    
    
    
    /// Implementation of \ref RunnerInterface::getRightPrecon()
    virtual PreconditionerInterface* getRightPrecon() override
    {
      return preconPair.r;
    }
    
  protected:
    LinearSolverInterface& oem;
    ITLSolver solver;
    
    /// create left/right preconditioners from parameters given in the init-file
    void setPrecon(PreconPair<MatrixType, VectorType>& pair)
    { FUNCNAME("ITL_Runner::setPrecon()");
    
      // Creator for the left preconditioner
      CreatorInterfaceName< ITL_PreconditionerBase<MatrixType, VectorType> >* leftCreator(NULL);

      // Creator for the right preconditioner
      CreatorInterfaceName< ITL_PreconditionerBase<MatrixType, VectorType> >* rightCreator(NULL);

      std::string preconType("no");
      std::string initFileStr = oem.getName() + "->left precon";
      Parameters::get(initFileStr, preconType);
      leftCreator = dynamic_cast<CreatorInterfaceName< ITL_PreconditionerBase<MatrixType, VectorType> >*>(
	CreatorMap<ITL_PreconditionerBase<MatrixType, VectorType> >::getCreator(preconType, initFileStr) );
      TEST_EXIT(leftCreator != NULL)
	("There is no creator for the given left preconditioner '%s'\n", preconType.c_str());
      leftCreator->setName(initFileStr);

      preconType = "no";
      initFileStr = oem.getName() + "->right precon";
      Parameters::get(initFileStr, preconType);
      rightCreator = dynamic_cast<CreatorInterfaceName< ITL_PreconditionerBase<MatrixType, VectorType> >*>(
	CreatorMap<ITL_PreconditionerBase<MatrixType, VectorType> >::getCreator(preconType, initFileStr) );
      TEST_EXIT(rightCreator != NULL)
	("There is no creator for the given right preconditioner '%s'\n", preconType.c_str());
      rightCreator->setName(initFileStr);

      pair.l = leftCreator->create();
      pair.r = rightCreator->create();
    }
    
  private:
    PreconPair<MatrixType, VectorType> preconPair;
  };

} // namespace AMDiS

#endif // AMDIS_ITL_RUNNER_H
