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



/** \file StandardProblemIteration.h */

#ifndef AMDIS_STANDARDPROBLEMITERATION_H
#define AMDIS_STANDARDPROBLEMITERATION_H

#include "AMDiS_fwd.h"
#include "ProblemIterationInterface.h"

namespace AMDiS {

  /// A master problem for a single non coupled problem.
  class StandardProblemIteration : public virtual ProblemIterationInterface
  {
  public:
    /// constructor
    StandardProblemIteration(ProblemStatBase *prob)
      : problem(prob)
    {}

    virtual ~StandardProblemIteration() {}

    /// Implementation of \ref ProblemIterationIterface::beginIteration()
    virtual void beginIteration(AdaptInfo *adaptInfo) override;

    /// Implementation of \ref ProblemIterationInterface::oneIteration()
    virtual Flag oneIteration(AdaptInfo *adaptInfo, Flag toDo) override;

    /// Implementation of \ref ProblemIterationInterface::endIteration()
    virtual void endIteration(AdaptInfo *adaptInfo) override;

    /// Implementation of \ref ProblemIterationInterface::getNumProblems()
    virtual int getNumProblems() override
    { 
      return 1; 
    }

    /// Implementation of \ref ProblemIterationInterface::getProblem(int)
    virtual ProblemStatBase *getProblem(int number = 0) override
    {
      return problem;
    }

    /// Returns the name of the problem. TODO: why is this virtual???
    virtual std::string getName() override;

    /// Nested assemblage and mesh adaption. TODO: make this prortected
    Flag buildAndAdapt(AdaptInfo *adaptInfo, Flag toDo);

  protected:
    /// The problem to solve.
    ProblemStatBase *problem;

    /// Info level
    static int info;
  };

}

#endif

