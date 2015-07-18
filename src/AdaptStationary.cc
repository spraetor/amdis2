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


#include "AdaptStationary.h"
#include "Initfile.h"
// #include "est/Estimator.h"
#include "ProblemIterationInterface.h"

#if HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/MeshDistributor.h"
#include "parallel/MpiHelper.h"
#endif

namespace AMDiS {

  AdaptStationary::AdaptStationary(std::string name,
				   ProblemIterationInterface& prob,
				   AdaptInfo& info) 
    : AdaptBase(name, &prob, &info)
  {
    Parameters::get(name + "->info", info);
  }


  int AdaptStationary::adapt()
  {
#if HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::MeshDistributor::globalMeshDistributor->initParallelization(); 
#endif

    // initial iteration
    if (adaptInfo->getSpaceIteration() == -1) {
      problemIteration->beginIteration(adaptInfo);
      problemIteration->oneIteration(adaptInfo, NO_ADAPTION);
      problemIteration->endIteration(adaptInfo);
      adaptInfo->incSpaceIteration();
    }

    // adaption loop
    while (!adaptInfo->spaceToleranceReached() && 
	   (adaptInfo->getSpaceIteration() < adaptInfo->getMaxSpaceIteration() || 
	    adaptInfo->getMaxSpaceIteration() < 0) ) {

      problemIteration->beginIteration(adaptInfo);
      Flag adapted = problemIteration->oneIteration(adaptInfo, FULL_ITERATION);
      problemIteration->endIteration(adaptInfo);

      int isAdapted = static_cast<bool>(adapted);
#if HAVE_PARALLEL_DOMAIN_AMDIS
      Parallel::mpi::globalAdd(isAdapted);
#endif
      if (isAdapted == 0)
	break;
      
      adaptInfo->incSpaceIteration();
    }

    return 0;
  }

} // end namespace AMDiS
