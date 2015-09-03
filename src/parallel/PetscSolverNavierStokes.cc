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



#include "parallel/PetscSolverNavierStokes.h"
#include "parallel/PetscHelper.h"
#include "TransformDOF.h"

namespace AMDiS
{
  namespace Parallel
  {

    using namespace std;


    PetscErrorCode pcSchurShell(PC pc, Vec x, Vec y)
    {
      void* ctx;
      PCShellGetContext(pc, &ctx);
      NavierStokesSchurData* data = static_cast<NavierStokesSchurData*>(ctx);

      // Project out constant null space
      {
        int vecSize;
        VecGetSize(x, &vecSize);
        PetscScalar vecSum;
        VecSum(x, &vecSum);
        vecSum = vecSum / static_cast<PetscScalar>(-1.0 * vecSize);
        VecShift(x, vecSum);
      }

      KSPSolve(data->kspLaplace, x, y);
      MatMult(data->matConDif, y, x);
      KSPSolve(data->kspMass, x, y);

      PetscFunctionReturn(0);
    }


    PetscSolverNavierStokes::PetscSolverNavierStokes(string name)
      : PetscSolverGlobalMatrix(name, false),
        pressureComponent(-1),
        pressureNullSpace(true),
        useOldInitialGuess(false),
        velocitySolutionMode(0),
        massSolutionMode(0),
        laplaceSolutionMode(0),
        massMatrixSolver(NULL),
        laplaceMatrixSolver(NULL),
        nu(NULL),
        invTau(NULL),
        solution(NULL),
        phase(NULL)
    {
      Parameters::get(initFileStr + "->navierstokes->pressure component",
                      pressureComponent);
      TEST_EXIT(pressureComponent >= 0)
      ("For using PETSc stokes solver you must define a pressure component!\n");

      Parameters::get(initFileStr + "->navierstokes->pressure null space",
                      pressureNullSpace);
      TEST_EXIT(pressureNullSpace)("This is not yet tested, may be wrong!\n");

      Parameters::get(initFileStr + "->navierstokes->use old initial guess",
                      useOldInitialGuess);

      Parameters::get(initFileStr + "->navierstokes->velocity solver",
                      velocitySolutionMode);

      Parameters::get(initFileStr + "->navierstokes->mass solver",
                      massSolutionMode);

      Parameters::get(initFileStr + "->navierstokes->laplace solver",
                      laplaceSolutionMode);
    }


    void PetscSolverNavierStokes::solvePetscMatrix(SystemVector& vec,
        AdaptInfo* adaptInfo)
    {
      FUNCNAME("PetscSolverNavierStokes::solvePetscMatrix()");

      if (useOldInitialGuess)
      {
        VecSet(getVecSolInterior(), 0.0);

        for (int i = 0; i < solution->getSize(); i++)
          setDofVector(getVecSolInterior(), solution->getDOFVector(i), i, true);

        vecSolAssembly();
        KSPSetInitialGuessNonzero(kspInterior, PETSC_TRUE);
      }

      PetscSolverGlobalMatrix::solvePetscMatrix(vec, adaptInfo);
    }


    void PetscSolverNavierStokes::initSolver(KSP& ksp)
    {
      FUNCNAME("PetscSolverNavierStokes::initSolver()");
      // Create FGMRES based outer solver

      MSG("CREATE POS 1: %p\n", &ksp);
      KSPCreate(domainComm, &ksp);
#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(ksp, getMatInterior(), getMatInterior());
#else
      KSPSetOperators(ksp, getMatInterior(), getMatInterior(), SAME_NONZERO_PATTERN);
#endif
      if (getInfo() >= 10)
        KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL);
      else if (getInfo() >= 20)
        KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, PETSC_NULL, PETSC_NULL);
      petsc_helper::setSolver(ksp, "ns_", KSPFGMRES, PCNONE, getRelative(), getTolerance(), getMaxIterations());

      // Create null space information.
      if (pressureNullSpace)
        setConstantNullSpace(ksp, pressureComponent, true);
    }


    void PetscSolverNavierStokes::initPreconditioner(PC pc)
    {
      FUNCNAME("PetscSolverNavierStokes::initPreconditioner()");

      Timer t;

      TEST_EXIT(nu)("nu pointer not set!\n");
      TEST_EXIT(invTau)("invtau pointer not set!\n");
      TEST_EXIT(solution)("solution pointer not set!\n");

      int dim = componentSpaces[pressureComponent]->getMesh()->getDim();

      vector<int> velocityComponents;
      for (int i = 0; i < dim; i++)
        velocityComponents.push_back(i);

      PCSetType(pc, PCFIELDSPLIT);
      PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
      PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_FULL);

      createFieldSplit(pc, "velocity", velocityComponents);
      createFieldSplit(pc, "pressure", pressureComponent);
      PCSetFromOptions(pc);

      KSPSetUp(kspInterior);

      KSP* subKsp;
      int nSubKsp;
      PCFieldSplitGetSubKSP(pc, &nSubKsp, &subKsp);


      TEST_EXIT(nSubKsp == 2)
      ("Wrong numer of KSPs inside of the fieldsplit preconditioner!\n");

      KSP kspVelocity = subKsp[0];
      KSP kspSchur = subKsp[1];
      PetscFree(subKsp);

      switch (velocitySolutionMode)
      {
      case 0:
        petsc_helper::setSolver(kspVelocity, "",
                                KSPRICHARDSON, PCHYPRE, 0.0, 1e-14, 1);
        break;
      case 1:
        petsc_helper::setSolverWithLu(kspVelocity, "", KSPPREONLY,
                                      PCLU, MATSOLVERMUMPS , 0.0, 1e-14, 1);
        break;
      default:
        ERROR_EXIT("No velocity solution mode %d available!\n", velocitySolutionMode);
      }


      KSPSetType(kspSchur, KSPPREONLY);
      PC pcSub;
      KSPGetPC(kspSchur, &pcSub);
      PCSetType(pcSub, PCSHELL);
      PCShellSetApply(pcSub, pcSchurShell);
      PCShellSetContext(pcSub, &matShellContext);

      if (pressureNullSpace)
        setConstantNullSpace(kspSchur);

      // === Mass matrix solver ===

      const FiniteElemSpace* pressureFeSpace = componentSpaces[pressureComponent];
      DOFMatrix massMatrix(pressureFeSpace, pressureFeSpace);
      {
        Operator massOp(pressureFeSpace, pressureFeSpace);
        ZeroOrderTerm* massTerm = NULL;
        if ((!phase) || (*nu == 0.0))
          massTerm = new Simple_ZOT;
        else
          massTerm = new VecAtQP_ZOT(phase);
        massOp.addTerm(massTerm);
        massMatrix.assembleOperator(massOp);
        delete massTerm;
      }
      massMatrixSolver = createSubSolver(pressureComponent, "mass_");
      massMatrixSolver->fillPetscMatrix(&massMatrix);


      // === Laplace matrix solver ===

      DOFMatrix laplaceMatrix(pressureFeSpace, pressureFeSpace);
      {
        Operator laplaceOp(pressureFeSpace, pressureFeSpace);
        SecondOrderTerm* laplaceTerm = NULL;
        if ((!phase) || (*nu == 0.0))
          laplaceTerm = new Simple_SOT;
        else
          laplaceTerm = new VecAtQP_SOT(phase);
        laplaceOp.addTerm(laplaceTerm);
        laplaceMatrix.assembleOperator(laplaceOp);
        delete laplaceTerm;
      }
      laplaceMatrixSolver = createSubSolver(pressureComponent, string("laplace_"));
      laplaceMatrixSolver->fillPetscMatrix(&laplaceMatrix);


      // === Create convection-diffusion operator ===

      DOFVector<double> vx(pressureFeSpace, "vx");
      DOFVector<double> vy(pressureFeSpace, "vy");
      DOFVector<double> vz(pressureFeSpace, "vz");
      DOFVector<double> vp(pressureFeSpace, "vp");
      vx.interpol(solution->getDOFVector(0));
      if (dim >= 2)
        vy.interpol(solution->getDOFVector(1));
      if (dim >= 3)
        vz.interpol(solution->getDOFVector(2));

      DOFMatrix conDifMatrix(pressureFeSpace, pressureFeSpace);

      {
        Operator conDifOp(pressureFeSpace, pressureFeSpace);

        ZeroOrderTerm* conDif0 = NULL;
        SecondOrderTerm* conDif1 = NULL;
        FirstOrderTerm* conDif2 = NULL, *conDif3 = NULL, *conDif4 = NULL;

        if (!phase)
        {
          MSG("INIT WITHOUT PHASE!\n");

          conDif0 = new Simple_ZOT(*invTau);
          conDifOp.addTerm(conDif0);
          conDif1 = new Simple_SOT(*nu);
          conDifOp.addTerm(conDif1);
          conDif2 = new VecAtQP_FOT(&vx, NULL, 0);
          conDifOp.addTerm(conDif2, GRD_PHI);
          if (dim >= 2)
          {
            conDif3 = new VecAtQP_FOT(&vy, NULL, 1);
            conDifOp.addTerm(conDif3, GRD_PHI);
          }
          if (dim == 3)
          {
            conDif4 = new VecAtQP_FOT(&vz, NULL, 2);
            conDifOp.addTerm(conDif4, GRD_PHI);
          }
        }
        else     // no phase given
        {

          vp.interpol(phase);

          if (*nu > 0.0)
          {
            conDif0 = new VecAtQP_ZOT(&vp, NULL, *invTau);
            conDifOp.addTerm(conDif0);
            conDif1 = new VecAtQP_SOT(&vp, NULL, *nu);
            conDifOp.addTerm(conDif1);
            conDif2 = new Vec2AtQP_FOT(&vx, &vp, NULL, 0);
            conDifOp.addTerm(conDif2, GRD_PHI);

            if (dim >= 2)
            {
              conDif3 = new Vec2AtQP_FOT(&vy, &vp, NULL, 1);
              conDifOp.addTerm(conDif3, GRD_PHI);
            }
            if (dim == 3)
            {
              conDif4 = new Vec2AtQP_FOT(&vz, &vp, NULL, 2);
              conDifOp.addTerm(conDif4, GRD_PHI);
            }
          }
          else
          {
            conDif0 = new VecAtQP_ZOT(&vp, new LinearInterpolation(*rho1,*rho2,*invTau));
            conDifOp.addTerm(conDif0);
            conDif1 = new VecAtQP_SOT(&vp, new LinearInterpolation(*nu1,*nu2));
            conDifOp.addTerm(conDif1);
            conDif2 = new Vec2AtQP_FOT(&vx, &vp, new LinearInterpolation2(*rho1,*rho2), 0);
            conDifOp.addTerm(conDif2, GRD_PHI);

            if (dim >= 2)
            {
              conDif3 = new Vec2AtQP_FOT(&vy, &vp, new LinearInterpolation2(*rho1,*rho2), 1);
              conDifOp.addTerm(conDif3, GRD_PHI);
            }
            if (dim == 3)
            {
              conDif4 = new Vec2AtQP_FOT(&vz, &vp, new LinearInterpolation2(*rho1,*rho2), 2);
              conDifOp.addTerm(conDif4, GRD_PHI);
            }
          }
        }/**/

        conDifMatrix.assembleOperator(conDifOp);

        delete conDif0;
        delete conDif1;
        delete conDif2;
        if (dim >= 2)
          delete conDif3;
        if (dim == 3)
          delete conDif4;
      }

      conDifMatrixSolver = createSubSolver(pressureComponent, "condif_");
      conDifMatrixSolver->fillPetscMatrix(&conDifMatrix);


      // === Setup solver ===

      matShellContext.kspMass = massMatrixSolver->getSolver();
      matShellContext.kspLaplace = laplaceMatrixSolver->getSolver();
      matShellContext.matConDif = conDifMatrixSolver->getMatInterior();

      switch (massSolutionMode)
      {
      case 0:
        petsc_helper::setSolver(matShellContext.kspMass, "mass_",
                                KSPCG, PCJACOBI, 0.0, 1e-14, 2);
        break;
      case 1:
        petsc_helper::setSolverWithLu(matShellContext.kspMass, "mass_", KSPRICHARDSON,
                                      PCLU, MATSOLVERMUMPS, 0.0, 1e-14, 1);
        break;
      default:
        ERROR_EXIT("No mass solution mode %d available!\n", massSolutionMode);
      }

      switch (laplaceSolutionMode)
      {
      case 0:
        petsc_helper::setSolver(matShellContext.kspLaplace, "laplace_",
                                KSPRICHARDSON, PCHYPRE, 0.0, 1e-14, 1);
        break;
      case 1:
        petsc_helper::setSolverWithLu(matShellContext.kspLaplace, "laplace_",
                                      KSPRICHARDSON, PCLU, MATSOLVERMUMPS,
                                      0.0, 1e-14, 1);
        break;
      default:
        ERROR_EXIT("No laplace solution mode %d available!\n", laplaceSolutionMode);
      }

      setConstantNullSpace(matShellContext.kspLaplace);

      MSG("Setup of Navier-Stokes preconditioner needed %.5f seconds\n",
          t.elapsed());
    }


    void PetscSolverNavierStokes::exitPreconditioner(PC pc)
    {
      FUNCNAME("PetscSolverNavierStokes::exitPreconditioner()");

      massMatrixSolver->destroyMatrixData();
      massMatrixSolver->destroyVectorData();
      laplaceMatrixSolver->destroyMatrixData();
      laplaceMatrixSolver->destroyVectorData();
      conDifMatrixSolver->destroyMatrixData();
      conDifMatrixSolver->destroyVectorData();

      massMatrixSolver->destroyVectorData();
      laplaceMatrixSolver->destroyVectorData();
      conDifMatrixSolver->destroyVectorData();


      delete massMatrixSolver;
      massMatrixSolver = NULL;

      delete laplaceMatrixSolver;
      laplaceMatrixSolver = NULL;

      delete conDifMatrixSolver;
      conDifMatrixSolver = NULL;
    }
  }
}
