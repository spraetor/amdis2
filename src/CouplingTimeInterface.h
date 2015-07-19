/** \file CouplingTimeInterface.h */

#pragma once

#include <vector>

#include "Flag.h"
#include "ProblemTimeInterface.h"

namespace AMDiS {

  /**
   * \ingroup Problem 
   *
   * \brief
   */
  class CouplingTimeInterface : public virtual ProblemTimeInterface
  {
  public:
    void addTimeInterface(ProblemTimeInterface *interface) 
    {
      problemInterfaces.push_back(interface);
    }

    /// Executes all needed operations when the simulation time changes.
    virtual void setTime(AdaptInfo *adaptInfo) override
    {
      for (ProblemTimeInterface* prob : problemInterfaces)
	       prob->setTime(adaptInfo);
    }

    /// Called at the beginning of each timestep
    virtual void initTimestep(AdaptInfo *adaptInfo) override 
    {
      for (ProblemTimeInterface* prob : problemInterfaces)
	       prob->initTimestep(adaptInfo);
    }

    /// Called at the end of each timestep.
    virtual void closeTimestep(AdaptInfo *adaptInfo) override 
    {
      for (ProblemTimeInterface* prob : problemInterfaces)
	       prob->closeTimestep(adaptInfo);
    }

    /// Solves the initial problem.
    virtual void solveInitialProblem(AdaptInfo *adaptInfo) override 
    {
      for (ProblemTimeInterface* prob : problemInterfaces)
	       prob->solveInitialProblem(adaptInfo);
    }

    /// Solves the initial problem.
    virtual void transferInitialSolution(AdaptInfo *adaptInfo) override 
    {
      for (ProblemTimeInterface* prob : problemInterfaces)
	       prob->transferInitialSolution(adaptInfo);
    }

  protected:
    /// vector of coupled time interfaces
    std::vector<ProblemTimeInterface*> problemInterfaces;
  };

} // end namespace AMDiS
