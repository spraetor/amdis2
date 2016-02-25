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



#include "parallel/PetscSolverCahnHilliard2.hpp"
#include "parallel/PetscHelper.hpp"
#include "TransformDOF.h"

namespace AMDiS
{
  namespace Parallel
  {


    using namespace std;


    PetscErrorCode pcChShell2b(PC pc, Vec b, Vec x) // solve Px=b
    {
      FUNCNAME("PCApply()");
      void* ctx;
      PCShellGetContext(pc, &ctx);
      CahnHilliardData2* data = static_cast<CahnHilliardData2*>(ctx);

      Vec y1, y2;
      VecDuplicate(b, &y1);
      VecDuplicate(b, &y2);

      KSPSolve(data->kspMplusK, b, y1);
      MatMult(data->matMass, y1, y2);
      KSPSolve(data->kspMplusK, y2, x);

      PetscFunctionReturn(0);
    }

    /// solve Cahn-Hilliard Preconditioner
    PetscErrorCode pcChSchurShell(PC pc, Vec b, Vec x) // solve Px=b
    {
      FUNCNAME("PCApply()");
      void* ctx;
      PCShellGetContext(pc, &ctx);
      CahnHilliardData2* data = static_cast<CahnHilliardData2*>(ctx);

      /// create S = M + eps^2*delta*K*D^(-1)*K where D=diag(M)
      Mat K, S;
      MatDuplicate(data->matMinusDeltaK, MAT_COPY_VALUES, &K);

      MatGetDiagonal(data->matMass, x);
      VecReciprocal(x);
      MatDiagonalScale(K, x, PETSC_NULL); 						// K' := -delta*D^(-1)*K
      MatMatMult(data->matMinusDeltaK, K, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &S); 		// S := -delta*K*K'
      MatAYPX(S, sqr(*data->eps)/(*data->delta), data->matMass, DIFFERENT_NONZERO_PATTERN); 	// S = eps^2/delta*S + M

      /// create new solver for S
      KSP kspS;
      KSPCreate(*data->mpiCommGlobal, &kspS);
#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(kspS, S, S);
#else
      KSPSetOperators(kspS, S, S, SAME_NONZERO_PATTERN);
#endif
      petsc_helper::setSolver(kspS, "S_", KSPFGMRES, PCSHELL, 1e-6, 1e-8, 1);
      {
        PC pc;
        KSPGetPC(kspS, &pc);
        PCShellSetApply(pc, pcChShell2b);
        PCShellSetContext(pc, data);
      }

      KSPSolve(kspS, b, x); 				// S*x2 = x1

      MatDestroy(&S);
      MatDestroy(&K);
      KSPDestroy(&kspS);

      PetscFunctionReturn(0);
    }



    PetscSolverCahnHilliard2::PetscSolverCahnHilliard2(string name)
      : PetscSolverGlobalMatrix(name),
        useOldInitialGuess(false),
        laplaceSolutionMode(0),
        massMatrixSolver(NULL),
        laplaceMatrixSolver(NULL),
        deltaKMatrixSolver(NULL),
        eps(NULL),
        delta(NULL),
        tau(NULL),
        solution(NULL),
        phase(NULL)
    {
      Parameters::get(initFileStr + "->use old initial guess",
                      useOldInitialGuess);
    }


    void PetscSolverCahnHilliard2::initSolver(KSP& ksp)
    {
      FUNCNAME("PetscSolverCahnHilliard2::initSolver()");

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
      petsc_helper::setSolver(ksp, "ch_", KSPFGMRES, PCNONE, getRelative(), getTolerance(), getMaxIterations());

      if (useOldInitialGuess)
        KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

    }


    void PetscSolverCahnHilliard2::initPreconditioner(PC pc)
    {
      FUNCNAME("PetscSolverCahnHilliard2::initPreconditioner()");

      MPI::COMM_WORLD.Barrier();
      double wtime = MPI::Wtime();


      if (tau)
      {
        delta = new double;
        *delta = sqrt(*tau);
      }

      vector<int> chPotentialComponent;
      chPotentialComponent.push_back(0);
      vector<int> chSchurComponent;
      chSchurComponent.push_back(1);

      PCSetType(pc, PCFIELDSPLIT);
      PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
      PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_FULL);

      createFieldSplit(pc, "ch_potential", chPotentialComponent);
      createFieldSplit(pc, "ch_schur", chSchurComponent);
      PCSetFromOptions(pc);

      KSPSetUp(kspInterior);

      KSP* subKspCH;
      int nSubKspCH;
      PCFieldSplitGetSubKSP(pc, &nSubKspCH, &subKspCH);


      TEST_EXIT(nSubKspCH == 2)
      ("Wrong numer of KSPs inside of the fieldsplit preconditioner!\n");

      KSP kspChPotential = subKspCH[0];
      KSP kspChSchur = subKspCH[1];
      PetscFree(subKspCH);

      KSPSetType(kspChSchur, KSPPREONLY);
      PC pcSub;
      KSPGetPC(kspChSchur, &pcSub);
      PCSetType(pcSub, PCSHELL);
      PCShellSetApply(pcSub, pcChSchurShell);
      PCShellSetContext(pcSub, &matShellContext);

      const FiniteElemSpace* feSpace = componentSpaces[0];

      // === Mass matrix solver ===
      DOFMatrix massMatrix(feSpace, feSpace);
      Operator massOp(feSpace, feSpace);
      Simple_ZOT zot;
      massOp.addTerm(&zot);
      massMatrix.assembleOperator(massOp);
      massMatrixSolver = createSubSolver(0, "mass_");
      massMatrixSolver->fillPetscMatrix(&massMatrix);

      // === Laplace matrix solver ===
      DOFMatrix laplaceMatrix(feSpace, feSpace);
      Operator laplaceOp(feSpace, feSpace);
      laplaceOp.addTerm(&zot); // M
      Simple_SOT sot2((*eps)*sqrt(*delta));
      laplaceOp.addTerm(&sot2); // eps*sqrt(delta)*K
      laplaceMatrix.assembleOperator(laplaceOp);
      laplaceMatrixSolver = createSubSolver(0, "MpK_");
      laplaceMatrixSolver->fillPetscMatrix(&laplaceMatrix);

      // === matrix (-delta*K) ===
      DOFMatrix deltaKMatrix(feSpace, feSpace);
      Operator laplaceOp2(feSpace, feSpace);
      Simple_SOT sot(-(*delta));
      laplaceOp2.addTerm(&sot); // -delta*K
      deltaKMatrix.assembleOperator(laplaceOp2);
      deltaKMatrixSolver = createSubSolver(0, "laplace_");
      deltaKMatrixSolver->fillPetscMatrix(&deltaKMatrix);


      // === Setup solver ===
      matShellContext.kspMplusK = laplaceMatrixSolver->getSolver();
      matShellContext.matMinusDeltaK = deltaKMatrixSolver->getMatInterior();
      matShellContext.eps = eps;
      matShellContext.delta = delta;
      matShellContext.kspMass = massMatrixSolver->getSolver();
      matShellContext.matMass = massMatrixSolver->getMatInterior();
      matShellContext.mpiCommGlobal = &(meshDistributor->getMpiComm(0));

      petsc_helper::setSolver(matShellContext.kspMass, "mass_", KSPCG, PCJACOBI, 0.0, 1e-14, 2);
      petsc_helper::setSolver(matShellContext.kspMplusK, "MpK_", KSPRICHARDSON, PCHYPRE, 0.0, 1e-14, 1);
      petsc_helper::setSolver(kspChPotential, "chPotential",  KSPRICHARDSON, PCHYPRE, 0.0, 1e-14, 1);

      MSG("Setup of Cahn-Hilliard preconditioner needed %.5f seconds\n",
          MPI::Wtime() - wtime);
    }


    void PetscSolverCahnHilliard2::exitPreconditioner(PC pc)
    {
      FUNCNAME("PetscSolverCahnHilliard2::exitPreconditioner()");

      massMatrixSolver->destroyMatrixData();
      massMatrixSolver->destroyVectorData();
      laplaceMatrixSolver->destroyMatrixData();
      laplaceMatrixSolver->destroyVectorData();
      deltaKMatrixSolver->destroyMatrixData();
      deltaKMatrixSolver->destroyVectorData();


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
