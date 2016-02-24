/** \file UmfPackSolver.h */

#pragma once

#if defined HAVE_UMFPACK && defined MTL_HAS_UMFPACK

#include <iostream>
#include <boost/numeric/mtl/operation/two_norm.hpp>
#include <boost/numeric/mtl/interface/umfpack_solve.hpp>
#include "solver/LinearSolver.h"

namespace AMDiS
{

  template <class MatrixType, class VectorType>
  struct Umfpack_Runner : public RunnerBase<MatrixType, VectorType>
  {
    using Super = RunnerBase<MatrixType, VectorType>;
    using Self  = Umfpack_Runner<MatrixType, VectorType>;

    /// Constructor.
    /**
     * Reads parameters for solver with name 'NAME':
     *   NAME->store symbolic       \ref store_symbolic
     *   NAME->symmetric strategy   \ref symmetric_strategy
     *   NAME->alloc init           \ref alloc_init
     **/
    Umfpack_Runner(LinearSolverInterface* oem_)
      : oem(*oem_),
        solver(NULL),
        store_symbolic(0),
        symmetric_strategy(0),
        alloc_init(0.7)
    {
      Parameters::get(oem.getName() + "->store symbolic", store_symbolic); // ?
      Parameters::get(oem.getName() + "->symmetric strategy", symmetric_strategy);
      Parameters::get(oem.getName() + "->alloc init", alloc_init);
    }

    /// Destructor.
    ~Umfpack_Runner()
    {
      if (solver != NULL)
      {
        delete solver;
        solver = NULL;
      }
    }

    /// Implementation of \ref RunnerBase::init()
    virtual void init(SolverMatrix<Matrix<DOFMatrix*>> const& /*A*/,
                      MatrixType const& fullMatrix) override
    {
      namespace umfpack = mtl::matrix::umfpack;
      if (solver != NULL)
      {
        delete solver;
        solver = NULL;
      }

      try
      {
        solver = new umfpack::solver<MatrixType>(fullMatrix, symmetric_strategy,
						 alloc_init);
      }
      catch (umfpack::error& e)
      {
        ERROR_EXIT("UMFPACK_ERROR(factorize, %d) = %s\n", e.code, e.what());
      }
    }

    /// Implementation of \ref RunnerBase::solve()
    virtual int solve(MatrixType const& A, VectorType& x,
		      VectorType const& b) override
    {
      FUNCNAME("Umfpack_Runner::solve()");
      TEST_EXIT(solver != NULL)("The umfpack solver was not initialized\n");
      namespace umfpack = mtl::matrix::umfpack;

      int code = 0;
      try
      {
        code = (*solver)(x, b);
      }
      catch (umfpack::error& e)
      {
        ERROR_EXIT("UMFPACK_ERROR(solve, %d) = %s\n", e.code, e.what());
      }

      VectorType r(b);
      r -= A * x;
      double residual = two_norm(r);
      oem.setResidual(residual);
      oem.setErrorCode(code);

      return code;
    }

    /// Implementation of \ref RunnerInterface::adjoint_solve()
    virtual void exit() override {}

  public:
    LinearSolverInterface& oem;

  private:
    mtl::matrix::umfpack::solver<MatrixType>* solver;

    int store_symbolic;

    int symmetric_strategy;

    double alloc_init;
  };


  /**
   * \ingroup Solver
   * \class AMDiS::UmfPackSolver
   * \brief \implements LinearSolver
   * Wrapper for the external UMFPACK solver:
   *   http://www.cise.ufl.edu/research/sparse/umfpack/
   *
   * This is a direct solver for large sparse matrices.
   */
  using UmfPackSolver = LinearSolver< MTLTypes::MTLMatrix, MTLTypes::MTLVector,
				      Umfpack_Runner< MTLTypes::MTLMatrix,
						      MTLTypes::MTLVector > >;

}

#endif // HAVE_UMFPACK
