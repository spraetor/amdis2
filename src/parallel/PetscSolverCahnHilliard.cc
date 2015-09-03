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



#include "PetscSolverCahnHilliard.h"
#include "parallel/PetscHelper.h"
#include "parallel/PetscSolverGlobalMatrix.h"

namespace AMDiS
{
  namespace Parallel
  {

    using namespace std;

    /// solve Cahn-Hilliard Preconditioner
    PetscErrorCode pcChShell(PC pc, Vec b, Vec x) // solve Px=b
    {
      void* ctx;
      PCShellGetContext(pc, &ctx);
      CahnHilliardData* data = static_cast<CahnHilliardData*>(ctx);

      Vec b1, b2, x1, x2;
      VecNestGetSubVec(b, 0, &b1);
      VecNestGetSubVec(b, 1, &b2);

      VecNestGetSubVec(x, 0, &x1);
      VecNestGetSubVec(x, 1, &x2);

      Vec y1, y2;
      VecDuplicate(b1, &y1);
      VecDuplicate(b2, &y2);

      //     MatGetDiagonal(data->matM, y2);
      //     VecReciprocal(y2);
      //     VecPointwiseMult(y1, y2, b1);
      KSPSolve(data->kspMass, b1, y1); 			// M*y1 = b1
      MatMultAdd(data->matMinusDeltaK, y1, b2, x1); 	// -> x1 := b2-delta*K*y1

      KSPSolve(data->kspLaplace, x1, y2); 			// (M+eps*sqrt(delta))*y2 = x1
      MatMult(data->matM, y2, x1); 			// x1 := M*y2

      KSPSolve(data->kspLaplace, x1, x2);			// (M+eps*sqrt(delta))*x2 = x1
      double factor = (*data->eps)/sqrt(*data->delta);
      VecCopy(x2, x1); 					// x1 := x2
      VecAXPBYPCZ(x1, 1.0, factor, -factor, y1, y2);	// x1 = 1*y1 + factor*y2 - factor*x1

      VecDestroy(&y1);
      VecDestroy(&y2);

      PetscFunctionReturn(0);
    }


    PetscSolverCahnHilliard::PetscSolverCahnHilliard(string name, double* epsPtr, double* deltaPtr)
      : PetscSolverGlobalBlockMatrix(name),
        massMatrixSolver(NULL),
        laplaceMatrixSolver(NULL),
        deltaKMatrixSolver(NULL),
        useOldInitialGuess(false),
        phase(NULL),
        eps(epsPtr),
        delta(deltaPtr),
        tau(NULL)
    {
      Parameters::get(initFileStr + "->use old initial guess", useOldInitialGuess);
    }

    void PetscSolverCahnHilliard::initSolver(KSP& ksp)
    {
      // Create FGMRES based outer solver
      KSPCreate(meshDistributor->getMpiComm(0), &ksp);
#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(ksp, getMatInterior(), getMatInterior());
#else
      KSPSetOperators(ksp, getMatInterior(), getMatInterior(), SAME_NONZERO_PATTERN);
#endif
      if (getInfo() >= 10)
        KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL);
      else if (getInfo() >= 20)
        KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, PETSC_NULL, PETSC_NULL);
      petsc_helper::setSolver(ksp, "ch_", KSPFGMRES, PCNONE, getRelative(), getTolerance(), getMaxIterations());
      KSPSetFromOptions(ksp);

      if (useOldInitialGuess)
        KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    }


    void PetscSolverCahnHilliard::initPreconditioner(PC pc)
    {
      FUNCNAME("PetscSolverCahnHilliard::initPreconditioner()");
      MSG("PetscSolverCahnHilliard::initPreconditioner()\n");

      if (tau)
      {
        delta = new double;
        *delta = (*gamma) * (*tau);
      }

      TEST_EXIT(eps && delta)("eps and/or delta pointers not set!\n");


      //     KSPSetUp(kspInterior);

      PCSetType(pc, PCSHELL);
      PCShellSetApply(pc, pcChShell);
      PCShellSetContext(pc, &matShellContext);

      const FiniteElemSpace* feSpace = componentSpaces[0];

      // === matrix M ===
      DOFMatrix laplaceMatrix(feSpace, feSpace);
      Operator laplaceOp(feSpace, feSpace);

      DOFMatrix massMatrix(feSpace, feSpace);
      Operator massOp(feSpace, feSpace);

      DOFMatrix deltaKMatrix(feSpace, feSpace);
      Operator laplaceOp2(feSpace, feSpace);

      if (phase)
      {
        VecAtQP_ZOT zot(phase, NULL);
        massOp.addTerm(&zot);
        laplaceOp.addTerm(&zot); // M
        VecAtQP_SOT sot2(phase, NULL, (*eps)*sqrt(*delta));
        laplaceOp.addTerm(&sot2); // eps*sqrt(delta)*K
        VecAtQP_SOT sot(phase, NULL, -(*delta));
        laplaceOp2.addTerm(&sot); // -delta*K
        massMatrix.assembleOperator(massOp);
        massMatrixSolver = createSubSolver(0, "mass_");
        massMatrixSolver->fillPetscMatrix(&massMatrix);

        // === matrix (M + eps*sqrt(delta)*K) ===
        laplaceMatrix.assembleOperator(laplaceOp);
        laplaceMatrixSolver = createSubSolver(0, "laplace_");
        laplaceMatrixSolver->fillPetscMatrix(&laplaceMatrix);

        // === matrix (-delta*K) ===
        deltaKMatrix.assembleOperator(laplaceOp2);
        deltaKMatrixSolver = createSubSolver(0, "laplace2_");
        deltaKMatrixSolver->fillPetscMatrix(&deltaKMatrix);


      }
      else
      {
        Simple_ZOT zot;
        massOp.addTerm(&zot);
        laplaceOp.addTerm(&zot); // M
        Simple_SOT sot2((*eps)*sqrt(*delta));
        laplaceOp.addTerm(&sot2); // eps*sqrt(delta)*K
        Simple_SOT sot(-(*delta));
        laplaceOp2.addTerm(&sot); // -delta*K

        massMatrix.assembleOperator(massOp);
        massMatrixSolver = createSubSolver(0, "mass_");
        massMatrixSolver->fillPetscMatrix(&massMatrix);

        // === matrix (M + eps*sqrt(delta)*K) ===
        laplaceMatrix.assembleOperator(laplaceOp);
        laplaceMatrixSolver = createSubSolver(0, "laplace_");
        laplaceMatrixSolver->fillPetscMatrix(&laplaceMatrix);

        // === matrix (-delta*K) ===
        deltaKMatrix.assembleOperator(laplaceOp2);
        deltaKMatrixSolver = createSubSolver(0, "laplace2_");
        deltaKMatrixSolver->fillPetscMatrix(&deltaKMatrix);
      }




      // === Setup solver ===
      matShellContext.kspMass = massMatrixSolver->getSolver();
      matShellContext.kspLaplace = laplaceMatrixSolver->getSolver();
      matShellContext.matM = massMatrixSolver->getMatInterior();
      matShellContext.matMinusDeltaK = deltaKMatrixSolver->getMatInterior();
      matShellContext.eps = eps;
      matShellContext.delta = delta;

      matShellContext.mpiCommGlobal= &(meshDistributor->getMpiComm(0));

      petsc_helper::setSolver(matShellContext.kspMass, "mass_", KSPCG, PCBJACOBI, 0.0, 1e-14, 2);
      //     petsc_helper::setSolver(matShellContext.kspMass, "mass_", KSPRICHARDSON, PCLU, 0.0, 1e-14, 1);
      //     {
      //             PC pc;
      //       KSPGetPC(matShellContext.kspMass, &pc);
      //       PCFactorSetMatSolverPackage(pc, MATSOLVERMUMPS);
      //     }

      petsc_helper::setSolver(matShellContext.kspLaplace, "laplace_", KSPRICHARDSON, PCHYPRE, 0.0, 1e-14, 2);
      //     petsc_helper::setSolver(matShellContext.kspLaplace, "laplace_", KSPRICHARDSON, PCLU, 0.0, 1e-14, 1);
      //     {
      //             PC pc;
      //       KSPGetPC(matShellContext.kspLaplace, &pc);
      //       PCFactorSetMatSolverPackage(pc, MATSOLVERMUMPS);
      //     }
      PCSetFromOptions(pc);
    }


    PetscSolver* PetscSolverCahnHilliard::createSubSolver(int component,
        string kspPrefix)
    {
      FUNCNAME("PetscSolverCahnHilliard::createSubSolver()");

      vector<const FiniteElemSpace*> fe;
      fe.push_back(componentSpaces[component]);

      PetscSolver* subSolver = new PetscSolverGlobalMatrix("");
      subSolver->setKspPrefix(kspPrefix.c_str());
      subSolver->setMeshDistributor(meshDistributor, 0);
      subSolver->init(fe, fe);

      ParallelDofMapping& subDofMap = subSolver->getDofMap();
      subDofMap[0] = dofMap[component];
      subDofMap.update();

      return subSolver;
    }

    void PetscSolverCahnHilliard::exitPreconditioner(PC pc)
    {
      FUNCNAME("PetscSolverNavierStokes::exitPreconditioner()");

      massMatrixSolver->destroyMatrixData();
      laplaceMatrixSolver->destroyMatrixData();
      deltaKMatrixSolver->destroyMatrixData();

      massMatrixSolver->destroyVectorData();
      laplaceMatrixSolver->destroyVectorData();
      deltaKMatrixSolver->destroyVectorData();


      delete massMatrixSolver;
      massMatrixSolver = NULL;

      delete laplaceMatrixSolver;
      laplaceMatrixSolver = NULL;

      delete deltaKMatrixSolver;
      deltaKMatrixSolver = NULL;

      if (tau)
      {
        delete delta;
      }
    }


  }
}
