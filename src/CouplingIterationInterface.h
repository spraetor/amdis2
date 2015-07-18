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



/** \file CouplingIterationInterface.h */

#ifndef AMDIS_COUPLINGITERATIONINTERFACE_H
#define AMDIS_COUPLINGITERATIONINTERFACE_H

#include "Flag.h"
#include "AdaptInfo.h"
#include "ProblemIterationInterface.h"

namespace AMDiS {

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
    virtual void beginIteration(AdaptInfo *adaptInfo);

    /** \brief
     * Determines the execution order of the single adaption steps. If adapt is
     * true, mesh adaption will be performed. This allows to avoid mesh adaption,
     * e.g. in timestep adaption loops of timestep adaptive strategies.
     */
    virtual Flag oneIteration(AdaptInfo *adaptInfo, Flag toDo = FULL_ITERATION);

    /// Called after each adaption loop iteration.
    virtual void endIteration(AdaptInfo *adaptInfo);

    virtual int getNumProblems();

    /// Returns number of managed problems
    virtual size_t getNumIterationInterfaces() { return problems.size(); }

    /** \brief
     * Returns the problem with the given number. If only one problem
     * is managed by this master problem, the number hasn't to be given.
     */
    virtual ProblemStatBase *getProblem(int number = 0);

    virtual ProblemIterationInterface *getIterationInterface(size_t number = 0) { return problems[number]; }

    /// Returns the name of the problem with the given number.
    virtual std::string getName(size_t number);
    virtual std::string getName() { return getName(0); }

    virtual void setSolveProblem(size_t number, bool flag = true) { solveProblem[number] = flag; }
    virtual void setSolveProblem(std::string name, bool flag = true);

    /// Function that serializes the problem plus information about the iteration.
    virtual void serialize(std::ostream &out) {};

    /// Function that deserializes the problem plus information about the iteration.
    virtual void deserialize(std::istream &in) {};

  protected:

    /// vector/map of coupled stationary problems
    std::vector<ProblemIterationInterface*> problems;
    std::vector<bool> solveProblem;
  };

}

#endif // AMDIS_COUPLINGITERATIONINTERFACE_H

