#pragma once

// AMDiS includes
#include "CreatorInterface.hpp"
#include "io/VtkWriter.hpp"
#include "nonlin/NonLinSolver.hpp"
#include "solver/LinearSolverInterface.hpp"

namespace AMDiS
{
  /**
   * \ingroup Solver
   *
   * \Brief
   * Implements the newton method for solving a non linear system. Sub class of
   * NonLinSolver.
   */
  class Newton : public NonLinSolver
  {
  public:
    /// Creator class used in the NonLinSolverMap.
    class Creator : public NonLinSolverCreator
    {
    public:
      virtual ~Creator() {}

      /// Returns a new Newton object.
      NonLinSolver* create()
      {
        return new Newton(this->name, this->linearSolver);
      }
    };

    /// Calls constructor of base class NonLinSolver
    Newton(const std::string& name, LinearSolverInterface* linSolver)
      : NonLinSolver(name, linSolver),
        b(NULL),
        buildCycle(1)
    {

      Parameters::get(name + "->build cycle", buildCycle);
    }

  private:
    /// Realisation of NonLinSolver::init
    void init() {}

    /// realisation of NonLinSolver::nlsolve
    int nlsolve(SolverMatrix<Matrix<DOFMatrix*>>& mat,
                SystemVector& x, SystemVector& rhs,
                AdaptInfo& adaptInfo,
                ProblemStat* prob)
    {
      FUNCNAME("Newton::nlsolve()");

      if (b == NULL)
        b = new SystemVector(x);

      double err = 0.0, errOld = -1.0;
      int iter, n;

      MSG("iter. |     this->residual |     red. |    n |\n");

      for (iter = 1; iter <= this->maxIter; iter++)
      {
        // Assemble DF(x) and F(x)
        if (iter == 1 || (buildCycle > 0 && (iter-1) % buildCycle == 0))
          prob->buildAfterCoarsen(adaptInfo, 0, true, true);
        else
          prob->buildAfterCoarsen(adaptInfo, 0, false, true);

        // Initial guess is zero
        b->set(0.0);

        // Solve linear system
        n = solveLinearSystem(mat, *b, rhs);

        // x = x + d
        x += *b;

        if (this->usedNorm == NO_NORM || this->usedNorm == L2_NORM)
          err = L2Norm(b);
        else
          err = H1Norm(b);


        if (iter == 1)
          this->initialResidual = err;

        if (errOld <= 0)
          MSG("%5d | %12.5e | -------- | %4d |\n", iter, err, n);
        else
          MSG("%5d | %12.5e | %8.2e | %4d |\n", iter, err, err/errOld, n);

        residual = err;
        if (err < this->tolerance)
        {
          MSG("Finished successfully!\n");
          return iter;
        }
        errOld = err;
      }

      MSG("iter. %d, residual: %12.5e\n", iter, err);
      MSG("tolerance %e not reached\n", this->tolerance);

      this->residual = err;

      return iter;
    }

    /// Realisation of NonLinSolver::exit
    void exit()
    {
      if (b != NULL)
      {
        delete b;
        b = NULL;
      }
    }

  private:
    /// Internal used data
    SystemVector* b;

    /// build matrix every ith iteration,
    /// 0...build matrix only once,
    /// i>=1...rebuild matrix in ith solver iteration,
    /// standard = 1
    int buildCycle;
  };

} // end namespace AMDiS
