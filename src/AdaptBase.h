/** \file AdaptBase.h */

#pragma once

#include <cstddef>
#include <string>

// #include "AMDiS_fwd.h"
// #include <traits/not_null.hpp>

namespace AMDiS
{
    // forward declarations
    class String;
    class ProblemIterationInterface;
    class ProblemTimeInterface;
    class AdaptInfo;


  /// Interface for adaption loops.
  class AdaptBase
  {
  public:
    /// Constructor
    AdaptBase(String const& sname,
            //   not_null<ProblemIterationInterface*> problemIteration_,
              ProblemIterationInterface* problemIteration_,
              AdaptInfo& adapt,
              ProblemTimeInterface* problemTime_ = NULL,
              AdaptInfo* initialAdaptInfo_ = NULL)
      : name(sname),
        problemIteration(problemIteration_),
        adaptInfo(adapt),
        problemTime(problemTime_),
        initialAdaptInfo(initialAdaptInfo_)
    { }

    /// Destructor
    virtual ~AdaptBase() {}

    /** \brief
     * Pure virtual method. Must be overloaded by sub classes to perform
     * a concrete adaption loop.
     */
    virtual int adapt() = 0;

    /// Returns \ref name
    String const& getName() const
    {
      return name;
    }

    /// Returns \ref problemIteration_
    ProblemIterationInterface* getProblemIteration() const
    {
      return problemIteration;
    }

    ///
    void setProblemIteration(ProblemIterationInterface* pii)
    {
      problemIteration = pii;
    }

    /// Returns \ref adaptInfo
    AdaptInfo& getAdaptInfo() const
    {
      return adaptInfo;
    }

    /// Returns \ref problemTime_
    ProblemTimeInterface* getProblemTime() const
    {
      return problemTime;
    }

    ///
    void setProblemTime(ProblemTimeInterface* pti)
    {
      problemTime = pti;
    }

    /// Returns \ref initialAdaptInfo_
    AdaptInfo& getInitialAdaptInfo() const
    {
      return *initialAdaptInfo;
    }

  protected:
    /// Name of the adaption loop
    String const& name;

    /// Problem iteration interface
    ProblemIterationInterface* problemIteration;

    /// Main adapt info
    AdaptInfo& adaptInfo;

    /// problem time interface
    ProblemTimeInterface* problemTime;

    /** \brief
     * Adapt info for initial adapt. Will be given to
     * problemTime_->solveInitialProblem().
     */
    AdaptInfo* initialAdaptInfo;

    /// Info level
    static int info;
  };

} // end namespace AMDiS
