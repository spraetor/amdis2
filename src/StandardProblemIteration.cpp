#include "StandardProblemIteration.hpp"

#include "AdaptInfo.hpp"
#include "Global.hpp"
#include "Initfile.hpp"
#include "ProblemStatBase.hpp"
#include "Timer.hpp"

namespace AMDiS
{
  int StandardProblemIteration::info = 10;


  ProblemStatBase& StandardProblemIteration::getProblem(int number)
  {
    FUNCNAME_DBG("StandardProblemIteration::getProblem");
    TEST_EXIT_DBG(number == 0)("Problem number out of range!\n");
    return problem;
  }

  ProblemStatBase& StandardProblemIteration::getProblem(std::string name)
  {
    FUNCNAME_DBG("StandardProblemIteration::getProblem");
    TEST_EXIT_DBG(name == problem.getName())("Problem name does not match!\n");
    return problem;
  }

  void StandardProblemIteration::beginIteration(AdaptInfo& adaptInfo)
  {
    FUNCNAME("StandardProblemIteration::beginIteration()");

    INFO(info, 4)("\n");
    INFO(info, 4)("begin of iteration number: %d\n", adaptInfo.getSpaceIteration() + 1);
    INFO(info, 4)("=============================\n");
  }


  Flag StandardProblemIteration::oneIteration(AdaptInfo& adaptInfo, Flag toDo)
  {
    Flag flag = buildAndAdapt(adaptInfo, toDo);

    if (toDo.isSet(SOLVE))
      problem.solve(adaptInfo, true, false);

    if (toDo.isSet(SOLVE_RHS))
      problem.solve(adaptInfo, true, false);

    if (toDo.isSet(ESTIMATE))
      problem.estimate(adaptInfo);

    return flag;
  }


  void StandardProblemIteration::endIteration(AdaptInfo& adaptInfo)
  {
    FUNCNAME("StandardProblemIteration::endIteration()");

    INFO(info, 4)("\n");
    INFO(info, 4)("end of iteration number: %d\n", adaptInfo.getSpaceIteration() + 1);
    INFO(info, 4)("=============================\n");
  }


  Flag StandardProblemIteration::buildAndAdapt(AdaptInfo& adaptInfo, Flag toDo)
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    FUNCNAME("StandardProblemIteration::buildAndAdapt()");
#endif

    Flag flag = 0, markFlag = 0;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Timer t;
#endif

    if (toDo.isSet(MARK))
      markFlag = problem.markElements(adaptInfo);
    else
      markFlag = 3;

    if (toDo.isSet(BUILD))
      problem.buildBeforeRefine(adaptInfo, markFlag);

    // refine
    if (toDo.isSet(ADAPT) && markFlag.isSet(MESH_REFINED))
      flag = problem.refineMesh(adaptInfo);

    if (toDo.isSet(BUILD))
      problem.buildBeforeCoarsen(adaptInfo, markFlag);

    // coarsen
    if (toDo.isSet(ADAPT) && markFlag.isSet(MESH_COARSENED))
      flag |= problem.coarsenMesh(adaptInfo);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    INFO(info, 8)("Local mesh adaption needed %.5f seconds\n", t.elapsed());
#endif

    if (toDo.isSet(BUILD))
      problem.buildAfterCoarsen(adaptInfo, markFlag, true, true);

    if (toDo.isSet(BUILD_RHS))
      problem.buildAfterCoarsen(adaptInfo, markFlag, false, true);

    return flag;
  }


  std::string StandardProblemIteration::getName() const
  {
    return problem.getName();
  }

} // end namespace AMDiS
