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



/** \file ProblemTimeInterface.h */

#ifndef AMDIS_PROBLEMTIMEINTERFACE_H
#define AMDIS_PROBLEMTIMEINTERFACE_H

#include "Flag.h"

namespace AMDiS {

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
    virtual void setTime(AdaptInfo *adaptInfo) = 0;

    /// Called at the beginning of each timestep
    virtual void initTimestep(AdaptInfo *adaptInfo) = 0;

    /// Called at the end of each timestep.
    virtual void closeTimestep(AdaptInfo *adaptInfo) = 0;

    /// Solves the initial problem.
    virtual void solveInitialProblem(AdaptInfo *adaptInfo) = 0;  

    /// Solves the initial problem.
    virtual void transferInitialSolution(AdaptInfo *adaptInfo) = 0;  

    /// Function that serializes the problem plus information about the iteration.
    virtual void serialize(std::ostream &out) = 0;

    /// Function that deserializes the problem plus information about the iteration.
    virtual void deserialize(std::istream &in) = 0;
  };

}

#endif
