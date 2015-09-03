/** \file CouplingIterationInterface.h */

#include "Flag.h"
#include "CouplingIterationInterface.h"

namespace AMDiS
{
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
  void CouplingIterationInterface::beginIteration(AdaptInfo* adaptInfo)
  {
    FUNCNAME("CouplingIterationInterface::beginIteration()");
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
  Flag CouplingIterationInterface::oneIteration(AdaptInfo* adaptInfo, Flag toDo)
  {
    Flag flag = 0;

    for (size_t i = 0; i < problems.size(); ++i)
    {
      if (solveProblem[i])
      {
        problems[i]->beginIteration(adaptInfo);
        flag |= problems[i]->oneIteration(adaptInfo, toDo);
        problems[i]->endIteration(adaptInfo);
      }
    }

    return flag;
  }


  /// Called after each adaption loop iteration.
  void CouplingIterationInterface::endIteration(AdaptInfo* adaptInfo)
  {
    FUNCNAME("CouplingIterationInterface::endIteration()");
    MSG("\n");
    MSG("end of iteration number: %d\n",
        adaptInfo->getTimestepNumber() + 1);
    MSG("==================================================\n");
  }


  /// Returns number of managed problems
  int CouplingIterationInterface::getNumProblems() const
  {
    int num= 0;
    for (size_t i = 0; i < problems.size(); ++i)
      num += problems[i]->getNumProblems();
    return num;
  }


  /** \brief
    * Returns the problem with the given number. If only one problem
    * is managed by this master problem, the number hasn't to be given.
    */
  ProblemStatBase* CouplingIterationInterface::getProblem(int number)
  {
    FUNCNAME("CouplingIterationInterface::getProblem");

    int maxNum = getNumProblems();
    TEST_EXIT(maxNum > number)("Problem number out of range.");

    int sum = 0;
    ProblemStatBase* probIter = NULL;
    for (size_t i = 0; i < problems.size(); ++i)
    {
      if (sum + problems[i]->getNumProblems() <= number)
        sum += problems[i]->getNumProblems();
      else
        probIter = problems[i]->getProblem(number - sum);
    }

    TEST_EXIT_DBG(probIter)("Problem not found. Should not happen, since number is in range.");
    return probIter;
  }


  /// Returns the name of the problem.
  std::string CouplingIterationInterface::getName(int number) const
  {
    TEST_EXIT(getNumIterationInterfaces() > number)("Problem number out of range.");

    return problems[number]->getName();
  }


  void CouplingIterationInterface::setSolveProblem(std::string name, bool flag)
  {
    for (size_t i = 0; i < problems.size(); ++i)
    {
      if (problems[i]->getName() == name)
      {
        setSolveProblem(i, flag);
        break;
      }
    }
  }

} // end namespace AMDiS
