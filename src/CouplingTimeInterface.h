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
      interfaces_.push_back(interface);
    }

    /// Executes all needed operations when the simulation time changes.
    virtual void setTime(AdaptInfo *adaptInfo) 
    {
      int size = static_cast<int>(interfaces_.size());
      for (int i = 0; i < size; i++)
	interfaces_[i]->setTime(adaptInfo);
    }

    /// Called at the beginning of each timestep
    virtual void initTimestep(AdaptInfo *adaptInfo) 
    {
      int size = static_cast<int>(interfaces_.size());
      for (int i = 0; i < size; i++)
	interfaces_[i]->initTimestep(adaptInfo);
    }

    /// Called at the end of each timestep.
    virtual void closeTimestep(AdaptInfo *adaptInfo) 
    {
      int size = static_cast<int>(interfaces_.size());
      for (int i = 0; i < size; i++)
	interfaces_[i]->closeTimestep(adaptInfo);
    }

    /// Solves the initial problem.
    virtual void solveInitialProblem(AdaptInfo *adaptInfo) 
    {
      int size = static_cast<int>(interfaces_.size());
      for (int i = 0; i < size; i++)
	interfaces_[i]->solveInitialProblem(adaptInfo);
    }

    /// Solves the initial problem.
    virtual void transferInitialSolution(AdaptInfo *adaptInfo) 
    {
      int size = static_cast<int>(interfaces_.size());
      for (int i = 0; i < size; i++)
	interfaces_[i]->transferInitialSolution(adaptInfo);
    }

    /// Function that serializes the problem plus information about the iteration.
    virtual void serialize(std::ostream &out) {};

    /// Function that deserializes the problem plus information about the iteration.
    virtual void deserialize(std::istream &in) {};

  protected:
    /// vector of coupled time interfaces
    std::vector<ProblemTimeInterface*> interfaces_;
  };

} // end namespace AMDiS
