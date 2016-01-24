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



/** \file PetscSolverNavierStokes.h */

#ifndef AMDIS_PETSC_SOLVER_NAVIER_STOKES_H
#define AMDIS_PETSC_SOLVER_NAVIER_STOKES_H

#include "parallel/PetscSolverGlobalMatrix.h"

namespace AMDiS
{
  namespace Parallel
  {

    struct NavierStokesSchurData
    {
      KSP kspMass;

      KSP kspLaplace;

      Mat matConDif;
    };

    class PetscSolverNavierStokes : public PetscSolverGlobalMatrix
    {
    private:

      class LinearInterpolation : public AbstractFunction<double, double>
      {
      public:
        LinearInterpolation(double c1, double c2, double factor=1.0)
          : AbstractFunction<double, double>(0)
        {
          a = (c1-c2)/2.0*factor;
          b = (c1+c2)/2.0*factor;
          cmin=std::min(c1,c2)*factor;
          cmax=std::max(c1,c2)*factor;
        }

        double operator()(const double& x) const
        {
          double result = b+a*x;
          if (result<cmin) result = cmin;
          if (result>cmax) result = cmax;
          return result;
        }
      private:
        double a,b,cmin,cmax;
      };


      class LinearInterpolation2 : public BinaryAbstractFunction<double, double, double>
      {
      public:
        LinearInterpolation2(double c1, double c2, double factor=1.0)
          : BinaryAbstractFunction<double, double, double>(0)
        {
          a = (c1-c2)/2.0*factor;
          b = (c1+c2)/2.0*factor;
          cmin=std::min(c1,c2)*factor;
          cmax=std::max(c1,c2)*factor;
        }

        double operator()(const double& u, const double& x) const
        {
          double result = b+a*x;
          if (result<cmin) result = cmin;
          if (result>cmax) result = cmax;
          return result * u;
        }
      private:
        double a,b,cmin,cmax;
      };

      class EinsMinus : public AbstractFunction<double, double>
      {
      public:
        EinsMinus(double d)
          : AbstractFunction<double, double>(0),
            c(d)
        {}

        double operator()(const double& x) const
        {
          return c * std::max(1.0-x,0.000001);
        }
      private:
        double c;
      };



    public:
      /// Creator class
      class Creator : public LinearSolverCreator
      {
      public:
        virtual ~Creator() {}

        /// Returns a new PetscSolver object.
        LinearSolverInterface* create()
        {
          return new PetscSolverNavierStokes(this->name);
        }
      };

      PetscSolverNavierStokes(std::string name);

      void solvePetscMatrix(SystemVector& vec, AdaptInfo& adaptInfo);

      void setStokesData(double* invTauPtr, SystemVector* vec, double* nu1_=NULL, double* nu2_=NULL, double* rho1_=NULL, double* rho2_=NULL)
      {
        nu = new double;
        (*nu) = 0.0;
        invTau = invTauPtr;
        solution = vec;
        nu1=nu1_;
        nu2=nu2_;
        rho1=rho1_;
        rho2=rho2_;
      }

      void setStokesData(double* nuPtr, double* invTauPtr, SystemVector* vec)
      {
        nu = nuPtr;
        invTau = invTauPtr;
        solution = vec;
      }


      void setPhase(DOFVector<double>* d)
      {
        phase = d;
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

      PetscSolver* massMatrixSolver, *laplaceMatrixSolver, *conDifMatrixSolver;

      NavierStokesSchurData matShellContext;

      double* nu, *invTau, *nu1,*nu2,*rho1,*rho2;

      SystemVector* solution;

      DOFVector<double>* phase;
    };

  }
}

#endif
