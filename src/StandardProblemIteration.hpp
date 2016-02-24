#pragma once

#include <string>

#include "AMDiS_fwd.hpp"
#include "ProblemIterationInterface.hpp"

namespace AMDiS
{
  /// A master problem for a single non coupled problem.
  class StandardProblemIteration : public virtual ProblemIterationInterface
  {
  public:
    /// constructor
    StandardProblemIteration(ProblemStatBase& prob)
      : problem(prob)
    {}

    /// Implementation of \ref ProblemIterationIterface::beginIteration()
    virtual void beginIteration(AdaptInfo& adaptInfo) override;

    /// Implementation of \ref ProblemIterationInterface::oneIteration()
    virtual Flag oneIteration(AdaptInfo& adaptInfo, Flag toDo) override;

    /// Implementation of \ref ProblemIterationInterface::endIteration()
    virtual void endIteration(AdaptInfo& adaptInfo) override;

    /// Returns the name of the problem.
    virtual std::string getName() const override;

    virtual int getNumProblems() const override
    {
      return 1;
    }

    /// Return the managed ProblemStat \ref problem, by number
    virtual ProblemStatBase& getProblem(int number = 0) override;

    /// Return the managed ProblemStat \ref problem, by name
    virtual ProblemStatBase& getProblem(std::string name) override;
    
  protected:
    /// Nested assemblage and mesh adaption.
    Flag buildAndAdapt(AdaptInfo& adaptInfo, Flag toDo);

  protected:
    /// The problem to solve.
    ProblemStatBase& problem;

    /// Info level
    static int info;
  };

} // end namespace AMDiS
