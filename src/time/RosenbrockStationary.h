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



/** \file RosenbrockStationary.h */

#ifndef AMDIS_ROSENBROCKSTATIONARY_H
#define AMDIS_ROSENBROCKSTATIONARY_H

#include "AMDiS_fwd.h"
#include "ProblemStat.h"
#include "SystemVector.h"
#include "time/RosenbrockMethod.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/ParallelProblemStat.h"
#endif

namespace AMDiS {

  struct RosenbrockBoundary 
  {
    RosenbrockBoundary(AbstractFunction<double, WorldVector<double> >* fct_,
		       AbstractFunction<double, WorldVector<double> >* fctDt_,
		       DOFVector<double> *vec_,
		       int row_, int col_)
      : fct(fct_), fctDt(fctDt_), vec(vec_), row(row_), col(col_) {}
    
    AbstractFunction<double, WorldVector<double> > *fct;
    AbstractFunction<double, WorldVector<double> > *fctDt; // dt[f](t,x)
    
    DOFVector<double> *vec;

    int row;
    int col;
  };

  /// Realization of a Rosenbrock time-discretization of M*d_t(X) = F[x]
  /** 
   * 1/(tau*gamma) M*Y_i^k  - J(t^k, X^k)[Y_i^k]  =  F[t_i^k, X_i^k] + sum_{j=1}^{i-1} c_ij / tau * M * Y_j^k
   *
   * with stageSolution[i]:
   * X_i^k = X^k + sum_{j=1}^{i-1} a_ij * Y_j^k
   *
   * oldTime: t^k, stageTime: t_i^k
   * 
   * and new solution
   * X^{k+1} = X^k + sum_{j=1}^{s} m_j * Y_j^k
   **/
  class RosenbrockStationary : public ProblemStat
  {
  public:
    RosenbrockStationary(std::string name, int componentShift_ = 0)
      : ProblemStat(name),
	first(true),
	componentShift(componentShift_),
	minusOne(-1.0),
	stageTime(0.0),
	oldTime(0.0),
	tauPtr(NULL),
	tauGamma(NULL),
	minusTauGamma(NULL),	
	invTauGamma(NULL),
	minusInvTauGamma(NULL)
    {}
    

    Flag oneIteration(AdaptInfo *adaptInfo, Flag toDo) override;
    
    virtual Flag stageIteration(AdaptInfo *adaptInfo, Flag flag, 
				bool asmMatrix, bool asmVector);
    
    virtual void estimateTimeError(AdaptInfo* adaptInfo);
    
    /// update solution vector and oldTime value
    void acceptTimestep(AdaptInfo* adaptInfo);
    
    /// Add operators of function F
    virtual void addOperator(Operator &op, int row, int col, 
			     double *factor = NULL, double *estFactor = NULL);

    /// Add operators of jacobian J = d_X(F)
    virtual void addJacobianOperator(Operator &op, int row, int col, 
				     double *factor = NULL, double *estFactor = NULL);

    virtual void addTimeOperator(int i, int j);

    // getting methods
    // _________________________________________________________________________
    
    DOFVector<double>* getUnVec(int i)
    {
      return unVec->getDOFVector(i);
    }

    DOFVector<double>* getStageSolution(int i)
    {
      return stageSolution->getDOFVector(i);
    }

    DOFVector<double>* getTimeRhsVec(int i)
    {
      return timeRhsVec->getDOFVector(i);
    }

    double* getTauGamma()
    {
      return tauGamma;
    }

    double* getTau()
    {
      return tauPtr;
    }

     double* getMinusTauGamma()
    {
      return minusTauGamma;
    }

    double* getInvTauGamma()
    {
      return invTauGamma;
    }

    double* getMinusInvTauGamma()
    {
      return minusInvTauGamma;
    }
    
    double* getStageTime()
    {
      return &stageTime;
    }
    
    double* getOldTime()
    {
      return &oldTime;
    }
    
    double* getTauGammaI()
    {
      return &tauGammaI;
    }

    // setting methods
    // _________________________________________________________________________
    
    void setRosenbrockMethod(RosenbrockMethod *method)
    {
      rm = method;
      init();
    }

    void setTau(double *ptr)
    {
      tauPtr = ptr;
    }

    void setTauGamma(double *ptr0, double *ptr1, double *ptr2, double *ptr3)
    {
      tauGamma = ptr0;
      minusTauGamma = ptr1;
      invTauGamma = ptr2;
      minusInvTauGamma = ptr3;
    }    
    
    void setOldTime(double t)
    {
      oldTime = t;
    }
    
    void setStageTime(double t)
    {
      stageTime = t;
    }


    // boundary conditions
    // _________________________________________________________________________

   /// Adds a Dirichlet boundary condition, where the rhs is given by an 
    /// abstract function.
    void addDirichletBC(BoundaryType type, int row, int col,
			AbstractFunction<double, WorldVector<double> > *fct) override;

    void addDirichletBC(BoundaryType type, int row, int col,
			AbstractFunction<double, WorldVector<double> > *fct,
			AbstractFunction<double, WorldVector<double> > *fctDt);
			
    /// Adds a Dirichlet boundary condition, where the rhs is given by a DOF
    /// vector.
    void addDirichletBC(BoundaryType type, int row, int col,
			DOFVector<double> *vec) override
    {
      FUNCNAME("RosenbrockStationary::addDirichletBC()");

      ERROR_EXIT("Not yet supported!\n");
    }

  protected:
    void init();

    void reset()
    {
      first = true;
      stageSolution->set(0.0);
      unVec->set(0.0);
      for (int i = 0; i < rm->getStages(); i++)
	stageSolutions[i]->set(0.0);
    }

  protected:    
    RosenbrockMethod *rm;

    SystemVector *stageSolution, *unVec, *timeRhsVec, *newUn, *tmp, *lowSol;

    std::vector<SystemVector*> stageSolutions;

    bool first;

    int componentShift;
    
    double minusOne;
    double stageTime; // t_n + alpha_i*tau
    double oldTime;
    double tauGammaI; // tau*gamma_i

    double *tauPtr;
    double *tauGamma, *minusTauGamma, *invTauGamma, *minusInvTauGamma;

    std::vector<RosenbrockBoundary> boundaries;
  };

}

#endif // AMDIS_ROSENBROCKSTATIONARY_H
