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


/** \file ParallelProblemStatBase.h */

#ifndef AMDIS_PARALLEL_PROBLEM_STAT_BASE_H
#define AMDIS_PARALLEL_PROBLEM_STAT_BASE_H

#include "parallel/MeshDistributor.hpp"
#include "ProblemStat.h"

namespace AMDiS
{
  namespace Parallel
  {

    /**
     * \ingroup Problem
     *
     * \brief
     * This class defines the stationary problem definition in parallel
     * computations.
     **/
    class ParallelProblemStat : public ProblemStatSeq
    {
    public:
      ParallelProblemStat(std::string nameStr,
                          ProblemIterationInterface* problemIteration = NULL);

      virtual ~ParallelProblemStat() {}

      void buildAfterCoarsen(AdaptInfo& adaptInfo, Flag flag,
                             bool assembleMatrix = true,
                             bool assembleVector = true);

      void initialize(Flag initFlag,
                      ProblemStatSeq* adoptProblem = NULL,
                      Flag adoptFlag = INIT_NOTHING);

      /// Must be called before Meshdistributor::initParallelization()
      /// is called.
      void addPeriodicBC(BoundaryType type, int row, int col);

    private:
      /// add PETSc and PMTL solvers to CreatorMap
      void addSolvers();

      /// add parallel PMTL preconditioners to CreatorMap
      void addPreconditioner();

    protected:
      MeshDistributor* meshDistributor;

    public:
      static double initTimeStamp;
      static bool initialized;
    };

  } // end namespace Parallel

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
  typedef detail::ProblemStat<Parallel::ParallelProblemStat> ProblemStat;
#endif

} // end namespace AMDiS

#endif
