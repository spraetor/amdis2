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
   * Implements the newton method with Armijo-rule for stepsize control
   * for solving a non linear system. Sub class of NonLinSolver.
   */
  class NewtonArmijo : public NonLinSolver
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
        return new NewtonArmijo(this->name, this->linearSolver);
      }
    };

    /// Calls constructor of base class NonLinSolver
    NewtonArmijo(const std::string& name, LinearSolverInterface* linSolver)
      : NonLinSolver(name, linSolver),
        b(NULL),
        buildCycle(1),
        delta(1.e-2),   // Abstiegsregulator, z.B. 1.e-2, 1.e-4
        alpha(0.5),     // Daempfungsfaktor, z.B. 0.5
        nLineSearch(5)  // maximale Anzahl line-seach Schritte
    {

      Parameters::get(name + "->build cycle", buildCycle);
      Parameters::get(name + "->armijo->delta", delta);
      Parameters::get(name + "->armijo->alpha", alpha);
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

      double err = 0.0, errOld = -1.0, lambda = 1.0;
      int iter, n;
      SystemVector x_test(x);

      MSG("iter. |     this->residual |     red. |    n |  lambda |\n");

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

        lambda = std::min(1.0, lambda/alpha);

        // x = x + d
        x_test = x + lambda * (*b);

        for (int k = 0; k < nLineSearch; ++k)
        {
          if (this->usedNorm == NO_NORM || this->usedNorm == L2_NORM)
            err = L2Norm(b);
          else
            err = H1Norm(b);
          // armijo rule
          if (err <= (1.0  - 2.0 * delta * lambda) * errOld)
            break;

          lambda *= alpha;
          x_test = x + lambda * (*b);
        }
        x = x_test;

        if (iter == 1)
          this->initialResidual = err;

        if (errOld <= 0)
          MSG("%5d | %12.5e | -------- | %4d | %5.2 |\n", iter, err, n, lambda);
        else
          MSG("%5d | %12.5e | %8.2e | %4d | %5.2 |\n", iter, err, err/errOld, n,lambda);

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

    /// Abstiegsregulator, z.B. 1.e-2, 1.e-4
    /// hinreichender Abstieg bei |f(x_{k+1]})|^2 < (1-2*delta*lambda_k)*|f(x_k})|^2
    double delta;

    /// Skalierungsfaktor, z.B. 0.5
    /// Anpassung des Daempfungsfaktors lambda um Faktor alpha:
    /// lambda_{k+1} = min(lambda_k/alpha, 1)
    double alpha;

    /// maximale Anzahl line-seach Schritte
    int nLineSearch;
  };

}
