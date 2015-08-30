#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

#include "ProblemInstat.h"
#include "AdaptStationary.h"
#include "AdaptInstationary.h"
#include "StandardProblemIteration.h"
#include "est/Estimator.h"
#include "io/FileWriter.h"

namespace AMDiS
{
  using namespace std;

  void ProblemInstatBase::solveInitialProblem(AdaptInfo* adaptInfo)
  {
    AdaptStationary initialAdapt((name + "->initial->adapt").c_str(),
                                 *(new StandardProblemIteration(initialProblem)),
                                 *adaptInfo);

    initialAdapt.adapt();
  }


  void ProblemInstat::transferInitialSolution(AdaptInfo* adaptInfo)
  {
    TEST_EXIT(adaptInfo->getTime() == adaptInfo->getStartTime())
    ("after initial solution: time != start time\n");
    problemStat->writeFiles(adaptInfo, true);
  }


  void ProblemInstat::closeTimestep(AdaptInfo* adaptInfo)
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    FUNCNAME("ProblemInstat::closeTimestep()");
#endif

    bool force = (adaptInfo->getTime() >= adaptInfo->getEndTime());
    problemStat->writeFiles(adaptInfo, force);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    MSG("Computational time for timestep: %.5f seconds\n",
        (MPI::Wtime() - lastTimepoint));
#endif
  }


  ProblemInstat::ProblemInstat(string sname,
                               ProblemStatSeq* prob,
                               ProblemStatBase* initialProb)
    : ProblemInstatBase(sname, initialProb),
      problemStat(prob),
      oldSolution(NULL)
  {}


  ProblemInstat::ProblemInstat(string sname,
                               ProblemStatSeq& prob)
    : ProblemInstatBase(sname, NULL),
      problemStat(&prob),
      oldSolution(NULL)
  {}


  ProblemInstat::ProblemInstat(string sname,
                               ProblemStatSeq& prob,
                               ProblemStatBase& initialProb)
    : ProblemInstatBase(sname, &initialProb),
      problemStat(&prob),
      oldSolution(NULL)
  {}


  ProblemInstat::~ProblemInstat()
  {
    if (oldSolution)
    {
      delete oldSolution;
      oldSolution = NULL;
    }
  }


  void ProblemInstat::initialize(Flag initFlag, ProblemInstat* adoptProblem,
                                 Flag adoptFlag)
  {
    FUNCNAME("ProblemInstat::initialize()");

    // === create vector for old solution ===

    if (oldSolution)
    {
      WARNING("oldSolution already created\n");
    }
    else
    {
      if (initFlag.isSet(INIT_UH_OLD))
        createUhOld();

      if (adoptProblem && adoptFlag.isSet(INIT_UH_OLD))
      {
        ProblemInstat* _adoptProblem = dynamic_cast<ProblemInstat*>(adoptProblem);
        TEST_EXIT(_adoptProblem)
        ("can't adopt oldSolution from problem which is not instationary and vectorial");
        TEST_EXIT(!oldSolution)("oldSolution already created");
        oldSolution = _adoptProblem->getOldSolution();
      }
    }

    if (!oldSolution)
      WARNING("no oldSolution created\n");
  }


  void ProblemInstat::createUhOld()
  {
    if (oldSolution)
    {
      WARNING("oldSolution already created\n");
    }
    else
    {
      int size = problemStat->getNumComponents();
      // create oldSolution
      oldSolution = new SystemVector("old solution", problemStat->getFeSpaces(), size);
      for (int i = 0; i < size; i++)
      {
        oldSolution->setDOFVector(i, new DOFVector<double>(problemStat->getFeSpace(i),
                                  name + "_uOld", true));
        oldSolution->getDOFVector(i)->setCoarsenOperation(COARSE_INTERPOL);

        if (problemStat->getEstimator(i))
          problemStat->getEstimator(i)->
          addUhOldToSystem(i, oldSolution->getDOFVector(i));
      }
    }
  }


  void ProblemInstat::initTimestep(AdaptInfo* adaptInfo)
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    lastTimepoint = MPI::Wtime();
#endif

    if (oldSolution)
      oldSolution->copy(*(problemStat->getSolution()));
  }

}
