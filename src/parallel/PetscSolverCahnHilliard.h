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



/** \file PetscSolverCahnHilliard.h */

#ifndef AMDIS_PETSC_SOLVER_CAHN_HILLIARD_H
#define AMDIS_PETSC_SOLVER_CAHN_HILLIARD_H

#include "parallel/PetscSolverGlobalBlockMatrix.h"

namespace AMDiS { namespace Parallel {

  struct CahnHilliardData {
    KSP kspMass, kspLaplace, kspLaplace2;
    Mat matM, matMinusDeltaK;
    double *eps, *delta; 
    MPI::Intracomm *mpiCommGlobal;

    
  };

  class PetscSolverCahnHilliard : public PetscSolverGlobalBlockMatrix
  {
  public:
    /// Creator class
    class Creator : public LinearSolverCreator
    {
    public:
      virtual ~Creator() {}

      /// Returns a new PetscSolver object.
      LinearSolverInterface* create() 
      { 
	return new PetscSolverCahnHilliard(this->name); 
      }
    };
    
    PetscSolverCahnHilliard(std::string name, double *epsPtr = NULL, double *deltaPtr = NULL);

    void setChData(double *epsPtr, double *deltaPtr, DOFVector<double> *p=NULL)
    {
      eps = epsPtr;
      delta = deltaPtr; // gamma*tau
      phase = p;
    }
    
    void setChData2(double *epsPtr, double *gammaPtr, double *tauPtr, DOFVector<double> *p=NULL)
    {
      eps = epsPtr;
      gamma = gammaPtr;
      tau = tauPtr;
      phase = p;
    }

  protected:
    void initSolver(KSP &ksp);

    void initPreconditioner(PC pc);
    void exitPreconditioner(PC pc);
    
    PetscSolver* createSubSolver(int component, std::string kspPrefix);

  private:
    int pressureComponent;

    PetscSolver *massMatrixSolver, *massMatrixSolver2, *laplaceMatrixSolver, *deltaKMatrixSolver;

    CahnHilliardData matShellContext;

    bool useOldInitialGuess;
    
    DOFVector<double> *phase;

    double *eps, *delta, *gamma, *tau;
  };

} }

#endif // AMDIS_PETSC_SOLVER_CAHN_HILLIARD_H

