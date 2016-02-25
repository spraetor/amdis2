#pragma once

// AMDiS includes
#include "AdaptInfo.hpp"
#include "Flag.hpp"
#include "ProblemIterationInterface.hpp"

namespace AMDiS
{

  /** \brief
   * Interface for master problems needed by the adaption loop. A master problem
   * can handle one single or multiple coupled problems. In the latter case,
   * the master problem can determine the execution order of the build, solve,
   * estimate, and adapt steps of the single problems in \ref oneIteration().
   * Standard execution order depends on the ordering of the problems in the
   * problem-list, that is filled by addProblem. Alternatively one can a access
   * each problem by an unique name.
   */
  class CouplingIterationInterface : public virtual ProblemIterationInterface
  {
  public:
    virtual ~CouplingIterationInterface() {}

    /// add problem by number
    virtual void addIterationInterface(ProblemIterationInterface* probIter, int number = -1);

    /// Called before each adaption loop iteration.
    virtual void beginIteration(AdaptInfo& adaptInfo) override;

    /** \brief
     * Determines the execution order of the single adaption steps. If adapt is
     * true, mesh adaption will be performed. This allows to avoid mesh adaption,
     * e.g. in timestep adaption loops of timestep adaptive strategies.
     */
    virtual Flag oneIteration(AdaptInfo& adaptInfo, Flag toDo = FULL_ITERATION) override;

    /// Called after each adaption loop iteration.
    virtual void endIteration(AdaptInfo& adaptInfo) override;

    virtual int getNumProblems() const override;

    /// Returns number of managed problems
    virtual int getNumIterationInterfaces() const
    {
      return static_cast<int>(problems.size());
    }

    /** \brief
     * Returns the problem with the given number. If only one problem
     * is managed by this master problem, the number hasn't to be given.
     */
    virtual ProblemStatBase& getProblem(int number = 0) override;

    virtual ProblemStatBase& getProblem(std::string name) override;

    virtual ProblemIterationInterface& getIterationInterface(int number = 0)
    {
      return *problems[number];
    }

    /// Returns the name of the problem with the given number.
    virtual std::string getName(int number) const;
    virtual std::string getName() const override
    {
      return "CouplingIterationInterface";
    }

    virtual void setSolveProblem(int number, bool flag = true)
    {
      solveProblem[number] = flag;
    }
    virtual void setSolveProblem(std::string name, bool flag = true);

  protected:
    /// vector/map of coupled stationary problems
    std::vector<ProblemIterationInterface*> problems;
    std::vector<bool> solveProblem;
  };

} // end namespace AMDiS
