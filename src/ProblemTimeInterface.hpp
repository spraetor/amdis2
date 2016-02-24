#pragma once

#include <iosfwd>

namespace AMDiS
{
  class AdaptInfo;

  // ============================================================================
  // ===== class ProblemTimeInterface ===========================================
  // ============================================================================

  /**
   * \ingroup Problem
   *
   * \brief
   * Interface for time dependent problems. Concrete problems must override
   * all pure virtual methods.
   */
  class ProblemTimeInterface
  {
  public:
    virtual ~ProblemTimeInterface() {};

    /// Called at the beginning of the adaption loop before any other call
    virtual void initTimeInterface() {};

    /// Executes all needed operations when the simulation time changes.
    virtual void setTime(AdaptInfo& adaptInfo) = 0;

    /// Called at the beginning of each timestep
    virtual void initTimestep(AdaptInfo& adaptInfo) = 0;

    /// Called at the end of each timestep.
    virtual void closeTimestep(AdaptInfo& adaptInfo) = 0;

    /// Solves the initial problem.
    virtual void solveInitialProblem(AdaptInfo& adaptInfo) = 0;

    /// Solves the initial problem.
    virtual void transferInitialSolution(AdaptInfo& adaptInfo) = 0;
  };

} // end namespace AMDiS
