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



/** \file NonLinSolver.h */

#ifndef AMDIS_NONLINSOLVER_H
#define AMDIS_NONLINSOLVER_H

#include <string>

#include "Global.h"
#include "CreatorInterface.h"
#include "MatrixVector.h"
#include "DOFMatrix.h"
#include "solver/LinearSolverInterface.h"
#include "ProblemStat.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/ParallelProblemStat.h"
#endif

namespace AMDiS
{

  /**
   * \ingroup Solver
   *
   * \brief
   * Pure virtual base class for specific non linear solvers.
   */
  class NonLinSolver
  {
  public:
    /** \brief
     * Constructor
     *
     * \param name     Name of this solver
     */
    NonLinSolver(const std::string& name, LinearSolverInterface* solver)
      : name(name),
        linSolver(solver),
        tolerance(1.e-8),
        maxIter(50),
        info(8),
        initialResidual(0.0),
        residual(0.0),
        usedNorm(NO_NORM)
    {
      Parameters::get(name + "->tolerance", tolerance);
      Parameters::get(name + "->max iteration", maxIter);
      Parameters::get(name + "->info", info);
      Parameters::get(name + "->norm", usedNorm);
    }

    /// Destructor
    virtual ~NonLinSolver() {}

    /** \brief
     * solves the non linear system. Uses sub class methods
     */
    inline int solve(SolverMatrix<Matrix<DOFMatrix*>>& mat,
                     SystemVector& solution, SystemVector& rhs,
                     AdaptInfo& adaptInfo,
                     ProblemStat* prob)
    {
      init();
      solution.set(0.0);
      int result = nlsolve(mat, solution, rhs, adaptInfo, prob);

      return result;
    }

    inline void setTolerance(double tol)
    {
      tolerance = tol;
    }

    inline double getTolerance()
    {
      return tolerance;
    }

    inline LinearSolverInterface* getLinearSolverInterface()
    {
      return linSolver;
    }

  protected:
    /// Allocates needed memory. Must be overriden in sub classes.
    virtual void init() = 0;

    /// Solves the non linear system. Must be overriden in sub classes.
    virtual int nlsolve(SolverMatrix<Matrix<DOFMatrix*>>& matVec,
                        SystemVector& x,
                        SystemVector& rhs,
                        AdaptInfo& adaptInfo,
                        ProblemStat* prob) = 0;

    /// Frees needed memory. Must be overriden in sub classes.
    virtual void exit() = 0;

    virtual int solveLinearSystem(SolverMatrix<Matrix<DOFMatrix*>>& mat,
                                  SystemVector& x,
                                  SystemVector& b)
    {
      FUNCNAME("NonLinSolver::solveLinearSystem()");
      TEST_EXIT(linSolver)("no solver\n");
      return linSolver->solveSystem(mat, x, b, true, false);
    }

  protected:
    /// Name of the solver.
    std::string name;

    /// Linear solver object.
    LinearSolverInterface* linSolver;

    /// Solver tolerance.
    double tolerance;

    /// Maximal number of iterations.
    int maxIter;

    /// Info level.
    int info;

    /// Initial residual.
    double initialResidual;

    /// Current residual.
    double residual;

    /// Used norm for convergent test.
    Norm usedNorm;
  };


  /// Interface for creators of concrete NonLinSolvers.
  class NonLinSolverCreator : public CreatorInterface<NonLinSolver>
  {
  public:
    virtual ~NonLinSolverCreator() {}

    void setName(std::string n)
    {
      name = n;
    }

    void setLinearSolverInterface(LinearSolverInterface* solver)
    {
      linearSolver = solver;
    }

  protected:
    std::string name;

    LinearSolverInterface* linearSolver;
  };

}

#include "Newton.h"
#include "NewtonArmijo.h"

#endif // AMDIS_NONLINSOLVER_H

