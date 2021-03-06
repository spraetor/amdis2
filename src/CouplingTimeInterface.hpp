#pragma once

// std c++ headers
#include <vector>

// AMDiS includes
#include "Flag.hpp"
#include "ProblemTimeInterface.hpp"

namespace AMDiS
{

  /**
   * \ingroup Problem
   *
   * \brief
   */
  class CouplingTimeInterface : public virtual ProblemTimeInterface
  {
  public:
    void addTimeInterface(ProblemTimeInterface* interface)
    {
      problemInterfaces.push_back(interface);
    }

    /// Executes all needed operations when the simulation time changes.
    virtual void setTime(AdaptInfo& adaptInfo) override
    {
      for (auto prob : problemInterfaces)
        prob->setTime(adaptInfo);
    }

    /// Called at the beginning of each timestep
    virtual void initTimestep(AdaptInfo& adaptInfo) override
    {
      for (auto prob : problemInterfaces)
        prob->initTimestep(adaptInfo);
    }

    /// Called at the end of each timestep.
    virtual void closeTimestep(AdaptInfo& adaptInfo) override
    {
      for (auto prob : problemInterfaces)
        prob->closeTimestep(adaptInfo);
    }

    /// Solves the initial problem.
    virtual void solveInitialProblem(AdaptInfo& adaptInfo) override
    {
      for (auto prob : problemInterfaces)
        prob->solveInitialProblem(adaptInfo);
    }

    /// Solves the initial problem.
    virtual void transferInitialSolution(AdaptInfo& adaptInfo) override
    {
      for (auto prob : problemInterfaces)
        prob->transferInitialSolution(adaptInfo);
    }

  protected:
    /// vector of coupled time interfaces
    std::vector<ProblemTimeInterface*> problemInterfaces;
  };

} // end namespace AMDiS
