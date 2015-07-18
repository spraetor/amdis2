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



/** \file ProblemIterationInterface.h */

#ifndef AMDIS_PROBLEMITERATIONINTERFACE_H
#define AMDIS_PROBLEMITERATIONINTERFACE_H

#include "Flag.h"

namespace AMDiS {

  class AdaptInfo;
  class ProblemStatBase;

  const Flag BUILD = 1;              // Assemble vectors and matrices
  const Flag BUILD_RHS = 2;          // Assemble rhs vectors only
  const Flag ADAPT = 4;              // Run adaption procedure
  const Flag SOLVE = 8;              // Solve system
  const Flag SOLVE_RHS = 16;         // Solve system, where only rhs vectors have changed
  const Flag ESTIMATE = 32;          // Estimate error
  const Flag MARK = 64;              // Mark elements

  const Flag FULL_ITERATION = BUILD | ADAPT | SOLVE | ESTIMATE | MARK;
  const Flag NO_ADAPTION = BUILD | SOLVE | ESTIMATE;

  /** \brief
   * Interface for master problems needed by the adaption loop. A master problem
   * can handle one single or multiple coupled problems. In the latter case,
   * the master problem can determine the execution order of the build, solve,
   * estimate, and adapt steps of the single problems in \ref oneIteration().
   */
  class ProblemIterationInterface
  {
  public:
    virtual ~ProblemIterationInterface() {}

    /// Called before each adaption loop iteration.
    virtual void beginIteration(AdaptInfo *adaptInfo) {}

    /** \brief
     * Determines the execution order of the single adaption steps. If adapt is
     * true, mesh adaption will be performed. This allows to avoid mesh adaption,
     * e.g. in timestep adaption loops of timestep adaptive strategies.
     */
    virtual Flag oneIteration(AdaptInfo *adaptInfo, Flag toDo = FULL_ITERATION) = 0;

    /// Called after each adaption loop iteration.
    virtual void endIteration(AdaptInfo *adaptInfo) {}

    /// Returns number of managed problems
    virtual int getNumProblems() = 0;

    /** \brief
     * Returns the problem with the given number. If only one problem
     * is managed by this master problem, the number hasn't to be given.
     */
    virtual ProblemStatBase *getProblem(int number = 0) = 0;

    /// Returns the problem with the given name. 
    virtual ProblemStatBase *getProblem(std::string name) 
    { 
      return NULL; 
    }

    /// Returns the name of the problem.
    virtual std::string getName() = 0;

    /// Function that serializes the problem plus information about the iteration.
    virtual void serialize(std::ostream &out) = 0;

    /// Function that deserializes the problem plus information about the iteration.
    virtual void deserialize(std::istream &in) = 0;
  };

}

#endif

