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


/** \file UmfPackSolver.h */

#ifndef AMDIS_UMFPACKSOLVER_H
#define AMDIS_UMFPACKSOLVER_H

#if defined HAVE_UMFPACK && defined MTL_HAS_UMFPACK

#include <iostream>
#include <boost/numeric/mtl/operation/two_norm.hpp>
#include <boost/numeric/mtl/interface/umfpack_solve.hpp>
#include "solver/LinearSolver.h"

namespace AMDiS
{

  template<typename MatrixType, typename VectorType>
  struct Umfpack_Runner : public RunnerBase<MatrixType, VectorType>
  {
    typedef RunnerBase<MatrixType, VectorType>      super;
    typedef Umfpack_Runner<MatrixType, VectorType>  self;

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

    /// Implementation of \ref RunnerBase::init()
    virtual void init(const SolverMatrix<Matrix<DOFMatrix*>>& A,
                      const MatrixType& fullMatrix) override
    {
      if (solver != NULL)
      {
        delete solver;
        solver = NULL;
      }

      try
      {
        solver = new mtl::matrix::umfpack::solver<MatrixType>(fullMatrix, symmetric_strategy, alloc_init);
      }
      catch (mtl::matrix::umfpack::error& e)
      {
        ERROR_EXIT("UMFPACK_ERROR(factorize, %d) = %s\n", e.code, e.what());
      }
    }

    /// Implementation of \ref RunnerBase::solve()
    virtual int solve(const MatrixType& A, VectorType& x, const VectorType& b) override
    {
      FUNCNAME("Umfpack_Runner::solve()");
      TEST_EXIT(solver != NULL)("The umfpack solver was not initialized\n");

      int code = 0;
      try
      {
        code = (*solver)(x, b);
      }
      catch (mtl::matrix::umfpack::error& e)
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

    ~Umfpack_Runner()
    {
      if (solver != NULL)
      {
        delete solver;
        solver = NULL;
      }
    }

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
  typedef LinearSolver<MTLTypes::MTLMatrix, MTLTypes::MTLVector, Umfpack_Runner<MTLTypes::MTLMatrix, MTLTypes::MTLVector>> UmfPackSolver;

}

#endif // HAVE_UMFPACK

#endif // AMDIS_UMFPACKSOLVER_H
