/** \file StandardProblemIteration.h */

#pragma once

#include "AMDiS_fwd.h"
#include "ProblemIterationInterface.h"

namespace AMDiS 
{
  /// A master problem for a single non coupled problem.
  class StandardProblemIteration : public virtual ProblemIterationInterface
  {
  public:
    /// constructor
    StandardProblemIteration(ProblemStatBase *prob)
      : problem(prob)
    {}

    virtual ~StandardProblemIteration() {}

    /// Implementation of \ref ProblemIterationIterface::beginIteration()
    virtual void beginIteration(AdaptInfo *adaptInfo) override;

    /// Implementation of \ref ProblemIterationInterface::oneIteration()
    virtual Flag oneIteration(AdaptInfo *adaptInfo, Flag toDo) override;

    /// Implementation of \ref ProblemIterationInterface::endIteration()
    virtual void endIteration(AdaptInfo *adaptInfo) override;

    /// Returns the name of the problem. TODO: why is this virtual???
    virtual std::string getName() override;
    
    virtual int getNumProblems() override
    {
      return 1;
    }
    
    virtual ProblemStatBase *getProblem(int number = 0) override
    { FUNCNAME_DBG("StandardProblemIteration::getProblem");
      TEST_EXIT_DBG(number == 0)("Problem number out of range!\n");
      return problem;
    }
    
    virtual ProblemStatBase *getProblem(std::string name) override
    { FUNCNAME_DBG("StandardProblemIteration::getProblem");
      TEST_EXIT_DBG(name == problem.getName())("Problem name does not match!\n");
      return problem;
    }

    /// Nested assemblage and mesh adaption. TODO: make this prortected
    Flag buildAndAdapt(AdaptInfo *adaptInfo, Flag toDo);

  protected:
    /// The problem to solve.
    ProblemStatBase *problem;

    /// Info level
    static int info;
  };

} // end namespace AMDiS
