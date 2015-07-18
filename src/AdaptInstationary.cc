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


#include "AdaptInstationary.h"
#include "Initfile.h"
#include "est/Estimator.h"
#include "ProblemIterationInterface.h"
#include "ProblemTimeInterface.h"
#include "Serializer.h"

#if HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/MeshDistributor.h"
#ifdef HAVE_PARALLEL_PETSC
#include <petsc.h>
#endif
#endif

using namespace std;

namespace AMDiS {

  AdaptInstationary::AdaptInstationary(std::string name,
				       ProblemIterationInterface &problemStat,  
				       AdaptInfo &info,
				       ProblemTimeInterface &problemInstat,
				       AdaptInfo &initialInfo,
				       time_t initialTimestampSet)
    : AdaptBase(name, &problemStat, &info, &problemInstat, &initialInfo),
      breakWhenStable(0)
  {
    strategy = 0;
    timeDelta1 = 0.7071;
    timeDelta2 = 1.4142;

    Parameters::get(name + "->strategy", strategy);
    Parameters::get(name + "->time delta 1", timeDelta1);
    Parameters::get(name + "->time delta 2", timeDelta2);
    Parameters::get(name + "->info", info);
    Parameters::get(name + "->break when stable", breakWhenStable);

    fixedTimestep = (info.getMinTimestep() == info.getMaxTimestep());
  }


  void AdaptInstationary::explicitTimeStrategy()
  {
    FUNCNAME("AdaptInstationary::explicitTimeStrategy()");

    // estimate before first adaption
    if (adaptInfo->getTime() <= adaptInfo->getStartTime())
      problemIteration->oneIteration(adaptInfo, ESTIMATE);


    // increment time
    adaptInfo->setTime(adaptInfo->getTime() + adaptInfo->getTimestep());

    problemTime->setTime(adaptInfo);

    INFO(info, 6)("time = %e, timestep = %e\n",
		  adaptInfo->getTime(), adaptInfo->getTimestep());

    adaptInfo->setSpaceIteration(0);
  
    // do the iteration
    problemIteration->beginIteration(adaptInfo);
    problemIteration->oneIteration(adaptInfo, FULL_ITERATION);
    problemIteration->endIteration(adaptInfo);
    adaptInfo->setLastProcessedTimestep(adaptInfo->getTimestep()); 
  }


  void AdaptInstationary::implicitTimeStrategy()
  {
    FUNCNAME("AdaptInstationary::implicitTimeStrategy()");

    do {
      adaptInfo->setTime(adaptInfo->getTime() + adaptInfo->getTimestep());
      problemTime->setTime(adaptInfo);

      INFO(info,6)("time = %e, try timestep = %e\n",
		   adaptInfo->getTime(), adaptInfo->getTimestep());

      problemIteration->oneIteration(adaptInfo, NO_ADAPTION);

      adaptInfo->incTimestepIteration();

      if (!fixedTimestep && 
	  !adaptInfo->timeToleranceReached() &&
	  adaptInfo->getTimestepIteration() <= adaptInfo->getMaxTimestepIteration() &&
	  !(adaptInfo->getTimestep() <= adaptInfo->getMinTimestep())) {
  
	adaptInfo->setTime(adaptInfo->getTime() - adaptInfo->getTimestep());
	adaptInfo->setTimestep(adaptInfo->getTimestep() * timeDelta1);
	continue;
      }


      adaptInfo->setSpaceIteration(0);


      // === Do space iterations only if the maximum is higher than 0. === 

      if (adaptInfo->getMaxSpaceIteration() > 0) {
    
	// === Space iterations. ===
	do {
	  problemIteration->beginIteration(adaptInfo);
	      
	  if (dbgMode) {
	    std::cout << "=== ADAPT INFO DEBUG MODE           ===\n";
	    std::cout << "=== in implicitTimeStrategy() ===\n";
	    std::cout << "=== space/time iteration "<< adaptInfo->getSpaceIteration()
	      <<" : "<< adaptInfo->getTimestepIteration() <<" ===\n";
	    adaptInfo->printTimeErrorLowInfo();
	  }
	  
	  Flag adapted = problemIteration->oneIteration(adaptInfo, FULL_ITERATION);
          int isAdapted = static_cast<bool>(adapted);
#if HAVE_PARALLEL_DOMAIN_AMDIS
	  AMDiS::Parallel::mpi::globalAdd(isAdapted)
#endif
	  if (isAdapted == 0) {
	    if (!fixedTimestep && 
		!adaptInfo->timeToleranceReached() &&
		!(adaptInfo->getTimestep() <= adaptInfo->getMinTimestep())) {
	      adaptInfo->setTime(adaptInfo->getTime() - adaptInfo->getTimestep());
	      adaptInfo->setTimestep(adaptInfo->getTimestep() * timeDelta2);
	      problemIteration->endIteration(adaptInfo);
	      adaptInfo->incSpaceIteration();
	      break;
	    }	
	  }

	  adaptInfo->incSpaceIteration();
	  problemIteration->endIteration(adaptInfo);
	  
	} while (!adaptInfo->spaceToleranceReached() && 
		 adaptInfo->getSpaceIteration() <= adaptInfo->getMaxSpaceIteration());

      } else {
	problemIteration->endIteration(adaptInfo);
      }


    } while(!adaptInfo->timeToleranceReached() &&
	    !(adaptInfo->getTimestep() <= adaptInfo->getMinTimestep()) && 
	    adaptInfo->getTimestepIteration() <= adaptInfo->getMaxTimestepIteration());  

    adaptInfo->setLastProcessedTimestep(adaptInfo->getTimestep()); 

    // After successful iteration/timestep the timestep will be changed according 
    // adaption rules for next timestep. 
    // First, check for increase of timestep
    if (!fixedTimestep && adaptInfo->timeErrorLow())
      adaptInfo->setTimestep(adaptInfo->getTimestep() * timeDelta2);

    // Second, check for decrease of timestep
    if (!fixedTimestep &&
	!adaptInfo->timeToleranceReached() &&
	!(adaptInfo->getTimestep() <= adaptInfo->getMinTimestep()))
	adaptInfo->setTimestep(adaptInfo->getTimestep() * timeDelta1);    
  }


  void AdaptInstationary::simpleAdaptiveTimeStrategy()
  {
    FUNCNAME("AdaptInstationary::simpleAdaptiveTimeStrategy()");

    // estimate before first adaption
    if (adaptInfo->getTime() <= adaptInfo->getStartTime())
      problemIteration->oneIteration(adaptInfo, ESTIMATE);

    adaptInfo->setTime(adaptInfo->getTime() + adaptInfo->getTimestep());
    problemTime->setTime(adaptInfo);
    
    INFO(info,6)("time = %e, timestep = %e\n",
		 adaptInfo->getTime(), adaptInfo->getTimestep());
    
    problemIteration->oneIteration(adaptInfo, FULL_ITERATION);

    adaptInfo->setLastProcessedTimestep(adaptInfo->getTimestep());
    
    // First, check for increase of timestep
    if (!fixedTimestep && adaptInfo->timeErrorLow())
      adaptInfo->setTimestep(adaptInfo->getTimestep() * timeDelta2);
    
    // Second, check for decrease of timestep
    if (!fixedTimestep &&
	!adaptInfo->timeToleranceReached() &&
	!(adaptInfo->getTimestep() <= adaptInfo->getMinTimestep()))
      adaptInfo->setTimestep(adaptInfo->getTimestep() * timeDelta1);    
  }


  void AdaptInstationary::oneTimestep()
  {
    FUNCNAME("AdaptInstationary::oneTimestep()");

    adaptInfo->setTimestepIteration(0);

    switch (strategy) {
    case 0:
      explicitTimeStrategy();
      break;
    case 1:
      implicitTimeStrategy();
      break;
    case 2:
      simpleAdaptiveTimeStrategy();
      break;
    default:
      ERROR_EXIT("Unknown strategy = %d!\n", strategy);
    }

    adaptInfo->incTimestepNumber();
  }


  int AdaptInstationary::adapt()
  {
    FUNCNAME("AdaptInstationary::adapt()");
    int errorCode = 0;

    TEST_EXIT(adaptInfo->getTimestep() >= adaptInfo->getMinTimestep())
      ("timestep < min timestep\n");
    TEST_EXIT(adaptInfo->getTimestep() <= adaptInfo->getMaxTimestep())
      ("timestep > max timestep\n");

    TEST_EXIT(adaptInfo->getTimestep() > 0)("timestep <= 0!\n");

#if HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::MeshDistributor::globalMeshDistributor->initParallelization(); 
#endif

    if (adaptInfo->getTimestepNumber() == 0) {
      adaptInfo->setTime(adaptInfo->getStartTime());
      initialAdaptInfo->setStartTime(adaptInfo->getStartTime());
      initialAdaptInfo->setTime(adaptInfo->getStartTime());

      problemTime->setTime(adaptInfo);

      // initial adaption
      problemTime->solveInitialProblem(initialAdaptInfo);
      problemTime->transferInitialSolution(adaptInfo);
    }

    while (!adaptInfo->reachedEndTime()) {
      iterationTimestamp = time(0);

      problemTime->initTimestep(adaptInfo);
      oneTimestep();
      problemTime->closeTimestep(adaptInfo);

      if (breakWhenStable && (adaptInfo->getSolverIterations() == 0)) 
	break;
    }

    return errorCode;
  }
}
