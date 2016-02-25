#include "nonlin/ProblemNonLin.hpp"

// AMDiS includes
#include "AdaptInfo.hpp"
#include "CreatorMap.hpp"
#include "Traverse.hpp"
#include "nonlin/NonLinSolver.hpp"

namespace AMDiS
{
  void ProblemNonLin::initialize(Flag initFlag,
                                 ProblemStatSeq* adoptProblem /*= NULL*/,
                                 Flag adoptFlag /*= INIT_NOTHING*/)
  {
    FUNCNAME("ProblemNonLin::initialize()");

    ProblemStat::initialize(initFlag, adoptProblem, adoptFlag);

    // === create nonlinear solver ===
    if (nonLinSolver)
    {
      WARNING("Nonlinear solver already created!\n");
    }
    else
    {
      if (initFlag.isSet(INIT_NONLIN_SOLVER))
        createNonLinSolver();

      if (adoptProblem && adoptFlag.isSet(INIT_NONLIN_SOLVER))
      {
        TEST_EXIT(nonLinSolver == NULL)("Nonlinear solver already created!\n");
        nonLinSolver = dynamic_cast<ProblemNonLin*>(adoptProblem)->getNonLinSolver();
      }
    }

    if (nonLinSolver == NULL)
      WARNING("No nonlinear solver created!\n");
  }


  void ProblemNonLin::createNonLinSolver()
  {
    // create non-linear solver
    std::string nonLinSolverType("no");
    std::string initFileStr(name + "->nonlin solver");
    Parameters::get(initFileStr, nonLinSolverType);

    NonLinSolverCreator* nonLinSolverCreator =
      dynamic_cast<NonLinSolverCreator*>(CreatorMap<NonLinSolver>::getCreator(nonLinSolverType, initFileStr));

    nonLinSolverCreator->setLinearSolverInterface(solver);
    nonLinSolverCreator->setName(name + "->nonlin solver");
    nonLinSolver = nonLinSolverCreator->create();
  }


  void ProblemNonLin::solve(AdaptInfo& adaptInfo, bool b0, bool b1)
  {
    TEST_EXIT(nonLinSolver)("no non-linear solver!\n");

    MSG("HERE A\n");

    nonLinSolver->solve(solverMatrix, *solution, *rhs, adaptInfo, this);
  }

}
