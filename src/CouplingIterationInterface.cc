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

#include "Flag.h"
#include "CouplingIterationInterface.h"

namespace AMDiS {

  /// add problem by number
  void CouplingIterationInterface::addIterationInterface(ProblemIterationInterface* probIter, int number)
  {
    size_t idx = std::min(problems.size(), static_cast<size_t>(number)); 
    if (number < 0)
      idx = problems.size();
    std::vector<ProblemIterationInterface*>::iterator pos = problems.begin() + idx;
    std::vector<bool>::iterator pos2 = solveProblem.begin() + idx;

    problems.insert(pos, probIter);
    solveProblem.insert(pos2, true);
  }


  /// Called before each adaption loop iteration.
  void CouplingIterationInterface::beginIteration(AdaptInfo *adaptInfo)
  { FUNCNAME("CouplingIterationInterface::beginIteration()");
    MSG("\n");
    int nTimesteps = (adaptInfo->getNumberOfTimesteps()
			? adaptInfo->getNumberOfTimesteps()
			: static_cast<int>((adaptInfo->getEndTime()-adaptInfo->getStartTime())/adaptInfo->getTimestep())
		     );
    MSG("begin of iteration number: %d/%d\n", 
      adaptInfo->getTimestepNumber() + 1,
      nTimesteps);
    MSG("==================================================\n");
  }


  /** \brief
    * Determines the execution order of the single adaption steps. If adapt is
    * true, mesh adaption will be performed. This allows to avoid mesh adaption,
    * e.g. in timestep adaption loops of timestep adaptive strategies.
    */
  Flag CouplingIterationInterface::oneIteration(AdaptInfo *adaptInfo, Flag toDo)
  {
    Flag flag = 0;

    for (size_t i = 0; i < problems.size(); ++i) {
      if (solveProblem[i]) {
        problems[i]->beginIteration(adaptInfo);
        flag |= problems[i]->oneIteration(adaptInfo, toDo);
        problems[i]->endIteration(adaptInfo);
      }
    }

    return flag;
  }


  /// Called after each adaption loop iteration.
  void CouplingIterationInterface::endIteration(AdaptInfo *adaptInfo)
  { FUNCNAME("CouplingIterationInterface::endIteration()");
    MSG("\n");
    MSG("end of iteration number: %d\n", 
      adaptInfo->getTimestepNumber() + 1);
    MSG("==================================================\n");
  }


  /// Returns number of managed problems
  int CouplingIterationInterface::getNumProblems()
  {
    size_t num= 0;
    for (size_t i = 0; i < problems.size(); ++i)
      num += problems[i]->getNumProblems();
    return num;
  }


  /** \brief
    * Returns the problem with the given number. If only one problem
    * is managed by this master problem, the number hasn't to be given.
    */
  ProblemStatBase *CouplingIterationInterface::getProblem(int number)
  {
    size_t maxNum = getNumProblems();
    if (maxNum <= static_cast<size_t>(number))
      throw(std::runtime_error("Problem number out of range."));
    
    size_t sum = 0;
    ProblemStatBase *probIter = NULL;
    for (size_t i = 0; i < problems.size(); ++i) {
      if (sum + problems[i]->getNumProblems() <= static_cast<size_t>(number))
        sum += problems[i]->getNumProblems();
      else
        probIter = problems[i]->getProblem(number - sum);
    }
    if (probIter == NULL)
      throw(std::runtime_error("Problem not found. Should not happen, since number is in range."));
    return probIter;
  }


  /// Returns the name of the problem.
  std::string CouplingIterationInterface::getName(size_t number)
  {
    if (static_cast<size_t>(problems.size()) <= number)
      throw(std::runtime_error("Problem number out of range."));

    return problems[number]->getName();
  }


  void CouplingIterationInterface::setSolveProblem(std::string name, bool flag)
  {
    for (size_t i = 0; i < problems.size(); ++i) {
      if (problems[i]->getName() == name) {
        setSolveProblem(i, flag);
        break;
      }
    }
  }

} // namespace AMDiS
