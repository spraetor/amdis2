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



/** \file ProblemInstat.h */

#ifndef AMDIS_PROBLEM_INSTAT_H
#define AMDIS_PROBLEM_INSTAT_H

#include "ProblemStat.h"
#include "ProblemTimeInterface.h"
#include "AdaptInstationary.h"

namespace AMDiS {

  /**
   * \ingroup Problem
   *
   * \brief
   * Base class for \ref ProblemInstat.
   */
  class ProblemInstatBase : public ProblemTimeInterface,
			    public ProblemStatBase
  {
  public:
    /// Constructor.
    ProblemInstatBase(std::string probName,
		      ProblemStatBase *initialProb)
      : name(probName),
	initialProblem(initialProb ? initialProb : this)
    {}

    /// Destructor.
    virtual ~ProblemInstatBase() {}

    /// Initialisation of the problem.
#if 0
    virtual void initialize(Flag initFlag,
			    ProblemInstat *adoptProblem = NULL,
			    Flag adoptFlag = INIT_NOTHING)
    {}
#endif

    virtual void setTime(AdaptInfo* adaptInfo) 
    {
      cTime = adaptInfo->getTime();
      tau = adaptInfo->getTimestep();
      invTau = 1.0 / tau;
    }

    virtual void solve(AdaptInfo* adaptInfo) {}

    virtual void solve(AdaptInfo *adaptInfo, bool, bool) 
    {
      solve(adaptInfo);
    }

    virtual void estimate(AdaptInfo *adaptInfo) {}

    virtual void buildBeforeRefine(AdaptInfo *adaptInfo, Flag) {}

    virtual void buildBeforeCoarsen(AdaptInfo *adaptInfo, Flag) {}

    virtual void buildAfterCoarsen(AdaptInfo *adaptInfo, Flag, bool, bool) {}

    virtual Flag markElements(AdaptInfo *adaptInfo) 
    { 
      return 0; 
    }

    virtual Flag refineMesh(AdaptInfo *adaptInfo) 
    { 
      return 0; 
    }

    virtual Flag coarsenMesh(AdaptInfo *adaptInfo) 
    { 
      return 0; 
    }

    /// Implementation of ProblemTimeInterface::closeTimestep().
    virtual void closeTimestep(AdaptInfo *adaptInfo) 
    {}

    /// Returns \ref name.
    inline std::string getName() 
    { 
      return name; 
    }

    /// Used by \ref problemInitial
    virtual void solveInitialProblem(AdaptInfo *adaptInfo);  

    double* getTime()
    {
      return &cTime;
    }

    double* getTau()
    {
      return &tau;
    }

    double* getInvTau()
    {
      return &invTau;
    }

    virtual void serialize(std::ostream &out) {}

    virtual void deserialize(std::istream &in) {}

  protected:
    /// Name of the problem.
    std::string name;

    ProblemStatBase *initialProblem;

    /// Time
    double cTime;

    /// Timestep
    double tau;

    /// 1 / timestep
    double invTau;
  };


  /**
   * \ingroup Problem
   *
   * \brief
   * Standard implementation of ProblemTimeInterface for a time
   * dependent problems.
   */
  class ProblemInstat : public ProblemInstatBase
  {
  public:
    /// Constructs a ProblemInstatVec with prob as its stationary problem.
    ProblemInstat(std::string name, 
		  ProblemStatSeq *prob,
		  ProblemStatBase *initialProb = NULL);

    ProblemInstat(std::string name, ProblemStatSeq &prob);

    ProblemInstat(std::string name, ProblemStatSeq &prob, ProblemStatBase &initialProb);

    /// Destructor.
    virtual ~ProblemInstat();

    /// Initialisation of the problem.
    virtual void initialize(Flag initFlag,
			    ProblemInstat *adoptProblem = NULL,
			    Flag adoptFlag = INIT_NOTHING);

    /// Used in \ref initialize().
    virtual void createUhOld();

    /// Implementation of ProblemTimeInterface::initTimestep().
    virtual void initTimestep(AdaptInfo *adaptInfo);

    /// Implementation of ProblemTimeInterface::closeTimestep().
    virtual void closeTimestep(AdaptInfo *adaptInfo);
  
    /// Returns \ref problemStat.
    inline ProblemStatSeq* getStatProblem() 
    {
      return problemStat; 
    }

    /// Returns \ref oldSolution.
    inline SystemVector *getOldSolution() 
    { 
      return oldSolution; 
    }

    inline DOFVector<double> *getOldSolution(int i)
    {
      return oldSolution->getDOFVector(i);
    }

    /// Used by \ref problemInitial
    virtual void transferInitialSolution(AdaptInfo *adaptInfo);  

  protected:
    /// Space problem solved in each timestep.
    ProblemStatSeq* problemStat;

    /// Solution of the last timestep.
    SystemVector *oldSolution;

    /// In parallel computations, we want to print the overall computational time
    /// that is used for one timestep.
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    double lastTimepoint;
#endif
  };

}

#endif
