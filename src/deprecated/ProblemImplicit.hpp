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



/** \file ProblemImplicit.h */

#ifndef AMDIS_PROBLEM_IMPLICIT_H
#define AMDIS_PROBLEM_IMPLICIT_H

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/ParallelProblemStat.hpp"
#else
#include "ProblemStat.h"
#endif
namespace AMDiS
{

  const Flag INIT_IMPLICIT_MESH = 0X2000L;

  class ProblemImplicit : public ProblemStat
  {
  public:
    ProblemImplicit(std::string name,
                    ProblemIterationInterface* problem = NULL)
      : ProblemStat(name, problem),
        r(0),
        phi1(0),
        phi2(0),
        levelSet(0)
    {}

    virtual ~ProblemImplicit();

    virtual void createMesh();

    virtual void initialize(Flag initFlag,
                            ProblemStatSeq* adoptProblem = NULL,
                            Flag adoptFlag = INIT_NOTHING);

    bool createImplicitMesh();

    DOFVector<double>* getSignedDistance(unsigned int im = 0,
                                         unsigned int m = 0) ;

    DOFVector<double>* getPhi1(unsigned int im = 0,
                               unsigned int m = 0) ;

    DOFVector<double>* getPhi2(unsigned int im = 0,
                               unsigned int m = 0) ;

    DOFVector<double>* getLevelset(unsigned int im = 0,
                                   unsigned int m = 0) ;

  protected:
    void readDofVec(const std::string&, DOFVector<double>*, Mesh*);

    void readR(const std::string&, double, Mesh*, int implMesh = 0, int comp = 0);

    void readPhi1(const std::string&, double, Mesh*, int implMesh = 0, int comp = 0);

    void readPhi2(const std::string&, double, Mesh*, int implMesh = 0, int comp = 0);

    std::string getDofFilename(std::string path, int implMesh = 0);

    double getEpsilon(std::string path, int implMesh = 0);

    int getType(std::string path, int implMesh = 0);

    bool createImplicitMesh(int p);

    bool createImplicitMesh(std::string, int, int);

    /// DOFVector for a signed distance
    std::vector<std::vector<DOFVector<double>*>> r;

    /// DOFVector for the phasefield function 0.5*(1-tanh(3*r/eps))
    std::vector<std::vector<DOFVector<double>*>> phi1;

    /// DOFVector for the phasefield function 0.5*(1+tanh(3*r/eps))
    std::vector<std::vector<DOFVector<double>*>> phi2;

    /// DOFVector for the levelset function
    /// (levelSet(x): x \in \Omega: 1, x \not \in Omega: -1, x \in \Gamma: 0)
    std::vector<std::vector<DOFVector<double>*>>  levelSet;

  };


}
#endif
