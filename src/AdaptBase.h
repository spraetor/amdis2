/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file AdaptBase.h */

#ifndef AMDIS_ADAPTBASE_H
#define AMDIS_ADAPTBASE_H

#include <string>
#include "AMDiS_fwd.h"
#include "Global.h"

namespace AMDiS {

  /// Interface for adaption loops.
  class AdaptBase
  {
  public:
    /// Constructor
    AdaptBase(std::string sname,
	      ProblemIterationInterface *problemIteration_,
	      AdaptInfo *adapt,
	      ProblemTimeInterface *problemTime_ = NULL,
	      AdaptInfo *initialAdaptInfo_ = NULL)
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
    std::string getName() const 
    { 
      return name; 
    }

    /// Returns \ref problemIteration_
    ProblemIterationInterface* getProblemIteration() const
    {
      return problemIteration;
    }

    ///
    void setProblemIteration(ProblemIterationInterface *pii) 
    {
      problemIteration = pii;
    }

    /// Returns \ref adaptInfo
    AdaptInfo* getAdaptInfo() const
    { 
      return adaptInfo; 
    }

    /// Returns \ref problemTime_
    ProblemTimeInterface* getProblemTime() const
    {
      return problemTime;
    }

    ///
    void setProblemTime(ProblemTimeInterface *pti) 
    {
      problemTime = pti;
    }

    /// Returns \ref initialAdaptInfo_
    AdaptInfo* getInitialAdaptInfo() const
    { 
      return initialAdaptInfo; 
    }

  protected:
    /// Name of the adaption loop
    std::string name;

    /// Problem iteration interface
    ProblemIterationInterface *problemIteration;

    /// Main adapt info
    AdaptInfo *adaptInfo;

    /// problem time interface
    ProblemTimeInterface *problemTime;

    /** \brief
     * Adapt info for initial adapt. Will be given to 
     * problemTime_->solveInitialProblem().
     */
    AdaptInfo *initialAdaptInfo;

    /// Info level
    static int info;
  };

}

#endif
