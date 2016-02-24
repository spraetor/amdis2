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



/** \file PetscSolverCahnHilliard2.h */

#ifndef AMDIS_PETSC_SOLVER_CAHN_HILLIARD2_H
#define AMDIS_PETSC_SOLVER_CAHN_HILLIARD2_H

#include "parallel/PetscSolverGlobalMatrix.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    struct CahnHilliardData2
    {
      KSP kspMass, kspMinusDeltaK, kspMplusK;
      Mat matMass, matMinusDeltaK;
      double* eps, *delta;
      MPI::Intracomm* mpiCommGlobal;
    };

    class PetscSolverCahnHilliard2 : public PetscSolverGlobalMatrix
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
          return new PetscSolverCahnHilliard2(this->name);
        }
      };

      PetscSolverCahnHilliard2(std::string name);

      void setChData(double* epsPtr, double* deltaPtr)
      {
        eps = epsPtr;
        delta = deltaPtr;
      }

      void setChData2(double* epsPtr, double* tauPtr, SystemVector* vec=NULL)
      {
        eps = epsPtr;
        tau = tauPtr;
        solution = vec;
      }


      void setPhase(DOFVector<double>* d, double eP3=0)
      {
        phase = d;
        epsPhase3 = eP3;
      }

    protected:
      void initSolver(KSP& ksp);

      void initPreconditioner(PC pc);

      void exitPreconditioner(PC pc);

    private:
      int pressureComponent;

      bool pressureNullSpace;

      /// If true, old solution is used for initial guess in solver phase.
      bool useOldInitialGuess;

      /// 0: approximate solve   1: direct solver
      int velocitySolutionMode;

      /// 0: approximate solve   1: direct solver
      int massSolutionMode;

      /// 0: approximate solve   1: direct solver
      int laplaceSolutionMode;

      PetscSolver* massMatrixSolver, *laplaceMatrixSolver, *deltaKMatrixSolver;

      CahnHilliardData2 matShellContext;

      double* eps, *delta, *tau;


      double epsPhase3;

      SystemVector* solution;

      DOFVector<double>* phase;

    };

  }
}

#endif
