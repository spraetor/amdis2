/** \file ITL_Runner.h */

#pragma once

#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/mtl/mtl.hpp>

#include "solver/LinearSolver.h"
#include "solver/ITL_Preconditioner.h"

namespace AMDiS
{
  template <class MatrixType, class VectorType>
  struct PreconPair
  {
    /// Pointer to the left preconditioner
    ITL_PreconditionerBase<MatrixType, VectorType>* l = NULL;

    /// Pointer to the right preconditioner
    ITL_PreconditionerBase<MatrixType, VectorType>* r = NULL;
    
    ~PreconPair()
    {
      if (l)
      {
        l->exit();
        delete l;
        l = NULL;
      }

      if (r)
      {
        r->exit();
        delete r;
        r = NULL;
      }
    }
  };


  /** \ingroup Solver
   *
   * \brief
   * Wrapper class for different MTL4 itl-solvers. These solvers
   * are parametrized by Matrix- and VectorType.
   **/
  template <class ITLSolver, class MatrixType, class VectorType>
  struct ITL_Runner : public RunnerBase<MatrixType, VectorType>
  {
    using Super = RunnerBase<MatrixType, VectorType>;

    /// Constructor.
    ITL_Runner(LinearSolverInterface* oemPtr)
      : oem(*oemPtr),
        solver(oem.getName())
    {
      setPrecon(preconPair);
    }

    /// Destructor.
    ~ITL_Runner() {}


    /// Implementation of \ref RunnerBase::init()
    virtual void init(SolverMatrix<Matrix<DOFMatrix*>> const& A,
                      MatrixType const& fullMatrix) override
    {
      preconPair.l->init(A, fullMatrix);
      preconPair.r->init(A, fullMatrix);
    }


    /// Implementation of \ref RunnerBase::solve()
    virtual int solve(MatrixType const& A , VectorType& x,
                      VectorType const& b) override
    {
      FUNCNAME("ITL_Runner::solve()");

      TEST_EXIT(preconPair.l)("there is no left preconditioner\n");
      TEST_EXIT(preconPair.r)("there is no right preconditioner\n");

      using value_type = typename mtl::Collection<MatrixType>::value_type;

      VectorType r(num_rows(b));
      r = b;
      if (two_norm(x) != 0)
      {
        r = A * x ;
        r -= b;
      }
      int error = 0;
      if (oem.getInfo() == 0)
      {
        // iteration that does not print residual per iteration
        itl::basic_iteration<value_type>
        iter(r, oem.getMaxIterations(), oem.getRelative(), oem.getTolerance());

        error = solver(A, x, b, *(preconPair.l), *(preconPair.r), iter);
        oem.setErrorCode(error);
        oem.setIterations(iter.iterations());
        oem.setResidual(iter.resid());
        oem.setRelativeResidual(iter.relresid());
      }
      else
      {
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
    virtual int adjoint_solve(MatrixType const& A ,
                              VectorType& x,
                              VectorType const& b) override
    {
      FUNCNAME("ITL_Runner::adjoint_solve()");
      ERROR_EXIT("Adjoint solve of preconditioner not yet implemented!\n");
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
    /// create left/right preconditioners from parameters given in the init-file
    void setPrecon(PreconPair<MatrixType, VectorType>& pair)
    {
      FUNCNAME("ITL_Runner::setPrecon()");

      // Creator for the left preconditioner
      CreatorInterfaceName<ITL_PreconditionerBase<MatrixType, VectorType>>* leftCreator(NULL);

      // Creator for the right preconditioner
      CreatorInterfaceName<ITL_PreconditionerBase<MatrixType, VectorType>>* rightCreator(NULL);

      std::string preconType("no");
      std::string initFileStr = oem.getName() + "->left precon";
      Parameters::get(initFileStr, preconType);
      leftCreator = dynamic_cast<CreatorInterfaceName<ITL_PreconditionerBase<MatrixType, VectorType>>*>(
                      CreatorMap<ITL_PreconditionerBase<MatrixType, VectorType>>::getCreator(preconType, initFileStr) );
      TEST_EXIT(leftCreator != NULL)
      ("There is no creator for the given left preconditioner '%s'\n", preconType.c_str());
      leftCreator->setName(initFileStr);

      preconType = "no";
      initFileStr = oem.getName() + "->right precon";
      Parameters::get(initFileStr, preconType);
      rightCreator = dynamic_cast<CreatorInterfaceName<ITL_PreconditionerBase<MatrixType, VectorType>>*>(
                       CreatorMap<ITL_PreconditionerBase<MatrixType, VectorType>>::getCreator(preconType, initFileStr) );
      TEST_EXIT(rightCreator != NULL)
      ("There is no creator for the given right preconditioner '%s'\n", preconType.c_str());
      rightCreator->setName(initFileStr);

      pair.l = leftCreator->create();
      pair.r = rightCreator->create();
    }

  protected:
    LinearSolverInterface& oem;
    ITLSolver solver;

  private:
    PreconPair<MatrixType, VectorType> preconPair;
  };

} // end namespace AMDiS
