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

#include "parallel/PetscSolverFetiOperators.hpp"
#include "parallel/PetscSolverFetiStructs.hpp"
#include "parallel/PetscSolverFetiTimings.hpp"

using namespace std;

namespace AMDiS
{
  namespace Parallel
  {

    void copyGlobalLocal(Vec globalB, Vec localB)
    {
      double* valueGlobal, *valueLocal;
      VecGetArray(globalB, &valueGlobal);
      VecGetArray(localB, &valueLocal);

      int nGlobal, nLocal;
      VecGetLocalSize(globalB, &nGlobal);
      VecGetLocalSize(localB, &nLocal);

      TEST_EXIT(nGlobal == nLocal)("Should not happen!\n");

      for (int i = 0; i < nGlobal; i++)
        valueLocal[i] = valueGlobal[i];

      VecRestoreArray(globalB, &valueGlobal);
      VecRestoreArray(localB, &valueLocal);
    }


    int petscMultMatSchurPrimal(Mat mat, Vec x, Vec y)
    {
      // S_PiPi = K_PiPi - K_PiB inv(K_BB) K_BPi

      void* ctx;
      MatShellGetContext(mat, &ctx);
      SchurPrimalData* data = static_cast<SchurPrimalData*>(ctx);

      MatMult(data->subSolver->getMatInteriorCoarse(), x, data->tmp_vec_b);
      data->subSolver->solveGlobal(data->tmp_vec_b, data->tmp_vec_b);
      MatMult(data->subSolver->getMatCoarseInterior(), data->tmp_vec_b,
              data->tmp_vec_primal);
      MatMult(data->subSolver->getMatCoarse(), x, y);
      VecAXPBY(y, -1.0, 1.0, data->tmp_vec_primal);

      return 0;
    }


    int petscMultMatSchurPrimalAugmented(Mat mat, Vec x, Vec y)
    {
      void* ctx;
      MatShellGetContext(mat, &ctx);
      SchurPrimalAugmentedData* data =
        static_cast<SchurPrimalAugmentedData*>(ctx);

      Vec x_primal, x_mu, y_primal, y_mu;

      if (data->nestedVec)
      {
        VecNestGetSubVec(x, 0, &x_primal);
        VecNestGetSubVec(x, 1, &x_mu);
        VecNestGetSubVec(y, 0, &y_primal);
        VecNestGetSubVec(y, 1, &y_mu);
      }
      else
      {
        VecDuplicate(data->tmp_vec_primal, &x_primal);
        VecDuplicate(data->tmp_vec_primal, &y_primal);

        PetscInt allLocalSize, allSize;
        VecGetLocalSize(x, &allLocalSize);
        VecGetSize(x, &allSize);

        PetscInt primalLocalSize, primalSize;
        VecGetLocalSize(x_primal, &primalLocalSize);
        VecGetSize(x_primal, &primalSize);

        TEST_EXIT_DBG(allSize > primalSize)("Should not happen!\n");
        TEST_EXIT_DBG(allLocalSize >= primalLocalSize)("Should not happen!\n");

        PetscInt muLocalSize = allLocalSize - primalLocalSize;
        PetscInt muSize = allSize - primalSize;

        VecCreateMPI(PETSC_COMM_WORLD, muLocalSize, muSize, &x_mu);
        VecCreateMPI(PETSC_COMM_WORLD, muLocalSize, muSize, &y_mu);

        PetscScalar* allValue;
        PetscScalar* primalValue;
        PetscScalar* muValue;
        VecGetArray(x, &allValue);
        VecGetArray(x_primal, &primalValue);
        VecGetArray(x_mu, &muValue);

        for (int i = 0; i < primalLocalSize; i++)
          primalValue[i] = allValue[i];
        for (int i = 0; i < muLocalSize; i++)
          muValue[i] = allValue[primalLocalSize + i];

        VecRestoreArray(x, &allValue);
        VecRestoreArray(x_primal, &primalValue);
        VecRestoreArray(x_mu, &muValue);
      }

      // inv(K_BB) K_BPi x_Pi
      MatMult(data->subSolver->getMatInteriorCoarse(), x_primal, data->tmp_vec_b0);
      data->subSolver->solveGlobal(data->tmp_vec_b0, data->tmp_vec_b0);

      // inv(K_BB) trans(J) trans(Q) x_mu
      MatMultTranspose(*(data->mat_augmented_lagrange), x_mu, data->tmp_vec_lagrange);
      MatMultTranspose(*(data->mat_lagrange), data->tmp_vec_lagrange, data->tmp_vec_b1);
      data->subSolver->solveGlobal(data->tmp_vec_b1, data->tmp_vec_b1);

      // y_Pi = (K_PiPi - K_PiB inv(K_BB) K_BPi) x_pi
      MatMult(data->subSolver->getMatCoarseInterior(), data->tmp_vec_b0,
              data->tmp_vec_primal);
      MatMult(data->subSolver->getMatCoarse(), x_primal, y_primal);
      VecAXPY(y_primal, -1.0, data->tmp_vec_primal);

      // y_Pi += (-K_PiB inv(K_BB) J^T Q^T) x_mu
      MatMult(data->subSolver->getMatCoarseInterior(), data->tmp_vec_b1,
              data->tmp_vec_primal);
      VecAXPY(y_primal, -1.0, data->tmp_vec_primal);

      // y_mu = (-Q J inv(K_BB) K_BPi) x_Pi + (-Q J inv(K_BB) J^T Q) x_mu
      //      = -Q J (inv(K_BB) K_BPi x_Pi + inv(K_BB) J^T Q x_mu)
      VecAXPY(data->tmp_vec_b0, 1.0, data->tmp_vec_b1);
      MatMult(*(data->mat_lagrange), data->tmp_vec_b0, data->tmp_vec_lagrange);
      MatMult(*(data->mat_augmented_lagrange), data->tmp_vec_lagrange, y_mu);
      VecScale(y_mu, -1.0);

      if (!data->nestedVec)
      {
        PetscInt allLocalSize;
        VecGetLocalSize(x, &allLocalSize);
        PetscInt primalLocalSize;
        VecGetLocalSize(x_primal, &primalLocalSize);
        PetscInt muLocalSize = allLocalSize - primalLocalSize;

        PetscScalar* allValue;
        PetscScalar* primalValue;
        PetscScalar* muValue;
        VecGetArray(y, &allValue);
        VecGetArray(y_primal, &primalValue);
        VecGetArray(y_mu, &muValue);

        for (int i = 0; i < primalLocalSize; i++)
          allValue[i] = primalValue[i];

        for (int i = 0; i < muLocalSize; i++)
          allValue[primalLocalSize + i] = muValue[i];

        VecRestoreArray(y, &allValue);
        VecRestoreArray(y_primal, &primalValue);
        VecRestoreArray(y_mu, &muValue);

        VecDestroy(&x_primal);
        VecDestroy(&y_primal);
        VecDestroy(&x_mu);
        VecDestroy(&y_mu);
      }

      return 0;
    }


    // y = mat * x
    int petscMultMatFeti(Mat mat, Vec x, Vec y)
    {
      FUNCNAME("petscMultMatFeti()");

      //    F = J inv(K_BB) trans(J) + J inv(K_BB) K_BPi inv(S_PiPi) K_PiB inv(K_BB) trans(J)
      // => F = J [I + inv(K_BB) K_BPi inv(S_PiPi) K_PiB] inv(K_BB) trans(J)

      double wtime = MPI::Wtime();

      void* ctx;
      MatShellGetContext(mat, &ctx);
      FetiData* data = static_cast<FetiData*>(ctx);

      MatMultTranspose(*(data->mat_lagrange), x, data->tmp_vec_b0);

      double wtime01 = MPI::Wtime();
      data->subSolver->solveGlobal(data->tmp_vec_b0, data->tmp_vec_b0);

      FetiTimings::fetiSolve01 += (MPI::Wtime() - wtime01);

      MatMult(*(data->mat_lagrange), data->tmp_vec_b0, data->tmp_vec_lagrange);

      MatMult(data->subSolver->getMatCoarseInterior(),
              data->tmp_vec_b0, data->tmp_vec_primal0);

      wtime01 = MPI::Wtime();
      KSPSolve(*(data->ksp_schur_primal), data->tmp_vec_primal0, data->tmp_vec_primal0);
      FetiTimings::fetiSolve02 += (MPI::Wtime() - wtime01);

      MatMult(data->subSolver->getMatInteriorCoarse(),
              data->tmp_vec_primal0, data->tmp_vec_b0);

      wtime01 = MPI::Wtime();
      data->subSolver->solveGlobal(data->tmp_vec_b0, data->tmp_vec_b0);
      FetiTimings::fetiSolve01 += (MPI::Wtime() - wtime01);

      MatMult(*(data->mat_lagrange), data->tmp_vec_b0, y);

      VecAXPBY(y, 1.0, 1.0, data->tmp_vec_lagrange);

      FetiTimings::fetiSolve += (MPI::Wtime() - wtime);

      return 0;
    }


    // y = mat * x
    int petscMultMatFetiInexact(Mat mat, Vec x, Vec y)
    {
      FUNCNAME("petscMultMatFetiInexact()");

      void* ctx;
      MatShellGetContext(mat, &ctx);
      FetiInexactData* data = static_cast<FetiInexactData*>(ctx);

      Vec xInterior, xPrimal, xLagrange;
      Vec yInterior, yPrimal, yLagrange;

      VecNestGetSubVec(x, 0, &xInterior);
      VecNestGetSubVec(x, 1, &xPrimal);
      VecNestGetSubVec(x, 2, &xLagrange);

      VecNestGetSubVec(y, 0, &yInterior);
      VecNestGetSubVec(y, 1, &yPrimal);
      VecNestGetSubVec(y, 2, &yLagrange);

      Vec tmpInterior;
      VecDuplicate(xInterior, &tmpInterior);
      MatMult(*(data->matBPi), xPrimal, tmpInterior);
      MatMultTranspose(*(data->mat_lagrange), xLagrange, yInterior);
      VecAXPY(yInterior, 1.0, tmpInterior);

      {
        {
          double* valueRhs, *valueTmp;
          VecGetArray(xInterior, &valueRhs);
          VecGetArray(data->tmp_vec_b0, &valueTmp);

          int nRhsLocal, nTmpLocal;
          VecGetLocalSize(xInterior, &nRhsLocal);
          VecGetLocalSize(data->tmp_vec_b0, &nTmpLocal);

          TEST_EXIT(nRhsLocal == nTmpLocal)("Should not happen!\n");

          for (int i = 0; i < nRhsLocal; i++)
            valueTmp[i] = valueRhs[i];

          VecRestoreArray(xInterior, &valueRhs);
          VecRestoreArray(data->tmp_vec_b0, &valueTmp);
        }

        MatMult(*(data->matBB), data->tmp_vec_b0, data->tmp_vec_b1);


        {
          double* valueRhs, *valueTmp;
          VecGetArray(tmpInterior, &valueRhs);
          VecGetArray(data->tmp_vec_b1, &valueTmp);

          int nRhsLocal, nTmpLocal;
          VecGetLocalSize(yInterior, &nRhsLocal);
          VecGetLocalSize(data->tmp_vec_b1, &nTmpLocal);

          TEST_EXIT(nRhsLocal == nTmpLocal)("Should not happen!\n");

          for (int i = 0; i < nRhsLocal; i++)
            valueRhs[i] = valueTmp[i];

          VecRestoreArray(tmpInterior, &valueRhs);
          VecRestoreArray(data->tmp_vec_b1, &valueTmp);
        }

        VecAXPY(yInterior, 1.0, tmpInterior);
      }

      {
        Vec tmpPrimal;
        VecDuplicate(xPrimal, &tmpPrimal);

        MatMult(*(data->matPiB), xInterior, tmpPrimal);
        MatMult(*(data->matPiPi), xPrimal, yPrimal);
        VecAXPY(yPrimal, 1.0, tmpPrimal);

        VecDestroy(&tmpPrimal);
      }

      MatMult(*(data->mat_lagrange), xInterior, yLagrange);


      VecDestroy(&tmpInterior);


      return 0;
    }


    PetscErrorCode pcInexactFetiShell(PC pc, Vec x, Vec y)
    {
      void* ctx;
      PCShellGetContext(pc, &ctx);
      FetiInexactPreconData* data = static_cast<FetiInexactPreconData*>(ctx);

      Vec xInterior, xPrimal, xLagrange;
      Vec yInterior, yPrimal, yLagrange;

      VecNestGetSubVec(x, 0, &xInterior);
      VecNestGetSubVec(x, 1, &xPrimal);
      VecNestGetSubVec(x, 2, &xLagrange);

      VecNestGetSubVec(y, 0, &yInterior);
      VecNestGetSubVec(y, 1, &yPrimal);
      VecNestGetSubVec(y, 2, &yLagrange);


      Vec tmpPrimal;
      VecDuplicate(xPrimal, &tmpPrimal);

      Vec tmpInterior0, tmpInterior1;
      VecDuplicate(xInterior, &tmpInterior0);
      VecDuplicate(xInterior, &tmpInterior1);

      {
        // === First Block ===

        copyGlobalLocal(xInterior, data->tmp_vec_b0);
        KSPSolve(data->ksp_interior, data->tmp_vec_b0, data->tmp_vec_b0);
        copyGlobalLocal(data->tmp_vec_b0, tmpInterior0);
        MatMult(*(data->matPiB), tmpInterior0, tmpPrimal);
        VecAYPX(tmpPrimal, -1.0, xPrimal);


        // === Second Block ===

        // tmpInterior0 already calculated

        KSPSolve(data->ksp_schur, tmpPrimal, tmpPrimal);


        // === Third Block ===

        // tmpPrimal already calculated

        MatMult(*(data->matBPi), tmpPrimal, tmpInterior1);
        copyGlobalLocal(tmpInterior1, data->tmp_vec_b0);
        KSPSolve(data->ksp_interior, data->tmp_vec_b0, data->tmp_vec_b0);
        copyGlobalLocal(data->tmp_vec_b0, tmpInterior1);

        VecAXPY(tmpInterior0, -1.0, tmpInterior1);

        VecCopy(tmpInterior0, yInterior);
        VecCopy(tmpPrimal, yPrimal);
      }


      {
        Vec tmpLagrange;
        VecDuplicate(xLagrange, &tmpLagrange);

        MatMult(*(data->mat_lagrange), yInterior, tmpLagrange);
        PCApply(data->pc_feti, tmpLagrange, yLagrange);

        PCApply(data->pc_feti, xLagrange, tmpLagrange);
        VecAXPY(yLagrange, -1.0, tmpLagrange);

        VecDestroy(&tmpLagrange);
      }

      VecDestroy(&tmpPrimal);
      VecDestroy(&tmpInterior0);
      VecDestroy(&tmpInterior1);

      PetscFunctionReturn(0);
    }


    // y = mat * x
    int petscMultMatFetiAugmented(Mat mat, Vec x, Vec y)
    {
      FUNCNAME("petscMultMatFetiAugmented()");

      void* ctx;
      MatShellGetContext(mat, &ctx);
      FetiData* data = static_cast<FetiData*>(ctx);

      Vec vec_mu0, vec_mu1;
      MatGetVecs(*(data->mat_augmented_lagrange), PETSC_NULL, &vec_mu0);
      VecDuplicate(vec_mu0, &vec_mu1);

      MatMultTranspose(*(data->mat_lagrange), x, data->tmp_vec_b0);
      data->subSolver->solveGlobal(data->tmp_vec_b0, data->tmp_vec_b0);

      MatMult(data->subSolver->getMatCoarseInterior(), data->tmp_vec_b0, data->tmp_vec_primal0);
      MatMult(*(data->mat_lagrange), data->tmp_vec_b0, data->tmp_vec_lagrange);
      MatMult(*(data->mat_augmented_lagrange), data->tmp_vec_lagrange, vec_mu0);

      Vec vec_array0[2] = {data->tmp_vec_primal0, vec_mu0};
      Vec vec_array1[2] = {data->tmp_vec_primal1, vec_mu1};
      Vec vec_nest0, vec_nest1;
      VecCreateNest(PETSC_COMM_WORLD, 2, PETSC_NULL, vec_array0, &vec_nest0);
      VecCreateNest(PETSC_COMM_WORLD, 2, PETSC_NULL, vec_array1, &vec_nest1);

      KSPSolve(*(data->ksp_schur_primal), vec_nest0, vec_nest1);

      // Step 1
      MatMult(*(data->mat_lagrange), data->tmp_vec_b0, y);

      // Step 2
      MatMult(data->subSolver->getMatInteriorCoarse(), data->tmp_vec_primal1, data->tmp_vec_b0);
      data->subSolver->solveGlobal(data->tmp_vec_b0, data->tmp_vec_b0);
      MatMult(*(data->mat_lagrange), data->tmp_vec_b0, data->tmp_vec_lagrange);
      VecAXPY(y, 1.0, data->tmp_vec_lagrange);

      // Step 3
      MatMultTranspose(*(data->mat_augmented_lagrange), vec_mu1, data->tmp_vec_lagrange);
      MatMultTranspose(*(data->mat_lagrange), data->tmp_vec_lagrange, data->tmp_vec_b0);
      data->subSolver->solveGlobal(data->tmp_vec_b0, data->tmp_vec_b0);
      MatMult(*(data->mat_lagrange), data->tmp_vec_b0, data->tmp_vec_lagrange);
      VecAXPY(y, 1.0, data->tmp_vec_lagrange);

      VecDestroy(&vec_mu0);
      VecDestroy(&vec_mu1);
      VecDestroy(&vec_nest0);
      VecDestroy(&vec_nest1);
      return 0;
    }


    int petscMultMatFetiInterface(Mat mat, Vec x, Vec y)
    {
      FUNCNAME("petscMultMatFetiInterface()");

      double wtime = MPI::Wtime();

      Vec x_interface, x_lagrange, y_interface, y_lagrange;
      VecNestGetSubVec(x, 0, &x_interface);
      VecNestGetSubVec(x, 1, &x_lagrange);
      VecNestGetSubVec(y, 0, &y_interface);
      VecNestGetSubVec(y, 1, &y_lagrange);

      void* ctx;
      MatShellGetContext(mat, &ctx);
      FetiData* data = static_cast<FetiData*>(ctx);


      // === Calculation of v_{B} ===

      // tmp_vec_b0 = J^{T} \lambda
      MatMultTranspose(*(data->mat_lagrange), x_lagrange, data->tmp_vec_b0);
      // tmp_vec_b1 = A_{B\Gamma} u_{\Gamma}
      MatMult(data->subSolver->getMatInteriorCoarse(1), x_interface, data->tmp_vec_b1);
      // tmp_vec_b0 = A_{B\Gamma} u_{\Gamma} + J^{T} \lambda
      VecAXPY(data->tmp_vec_b0, 1.0, data->tmp_vec_b1);


      // == Calculation of v_{\Pi}

      // tmp_vec_primal0 = A_{\Pi\Gamma} u_{\Gamma}
      MatMult(data->subSolver->getMatCoarse(0, 1), x_interface, data->tmp_vec_primal0);


      // === Calculate action of FETI-DP operator ===

      double wtime01 = MPI::Wtime();
      // tmp_vec_b0 = A_{BB}^{-1} v_{B}
      data->subSolver->solveGlobal(data->tmp_vec_b0, data->tmp_vec_b0);
      FetiTimings::fetiSolve01 += (MPI::Wtime() - wtime01);

      // tmp_vec_primal1 = A_{\Pi B} A_{BB}^{-1} v_{B}
      MatMult(data->subSolver->getMatCoarseInterior(),
              data->tmp_vec_b0, data->tmp_vec_primal1);

      // tmp_vec_primal0 = v_{\Pi} - A_{\Pi B} A_{BB}^{-1} v_{B}
      VecAXPY(data->tmp_vec_primal0, -1.0, data->tmp_vec_primal1);

      wtime01 = MPI::Wtime();
      // tmp_vec_primal0 = S_{\Pi\Pi}^{-1} (v_{\Pi} - A_{\Pi B} A_{BB}^{-1} v_{B})
      KSPSolve(*(data->ksp_schur_primal), data->tmp_vec_primal0, data->tmp_vec_primal0);
      FetiTimings::fetiSolve02 += (MPI::Wtime() - wtime01);

      // tmp_vec_b1 = A_{B\Pi} S_{\Pi\Pi}^{-1} (v_{\Pi} - A_{\Pi B} A_{BB}^{-1} v_{B})
      MatMult(data->subSolver->getMatInteriorCoarse(),
              data->tmp_vec_primal0, data->tmp_vec_b1);

      wtime01 = MPI::Wtime();
      // tmp_vec_b1 = A_{BB}^{-1} A_{B\Pi} S_{\Pi\Pi}^{-1} (v_{\Pi} - A_{\Pi B} A_{BB}^{-1} v_{B})
      data->subSolver->solveGlobal(data->tmp_vec_b1, data->tmp_vec_b1);
      FetiTimings::fetiSolve01 += (MPI::Wtime() - wtime01);

      // tmp_vec_b0 = A_{BB}^{-1} v_{B} - A_{BB}^{-1} A_{B\Pi} S_{\Pi\Pi}^{-1} (v_{\Pi} - A_{\Pi B} A_{BB}^{-1} v_{B})
      VecAXPY(data->tmp_vec_b0, -1.0, data->tmp_vec_b1);


      // === Calculate projection to interface and constraint variable ===

      // y_interface = A_{\Gamma B} tmp_vec_b0
      MatMult(data->subSolver->getMatCoarseInterior(1), data->tmp_vec_b0, y_interface);

      // tmp_vec_primal0 = S_{\Pi\Pi}^{-1} (v_{\Pi} - A_{\Pi B} A_{BB}^{-1} v_{B})
      // tmp_vec_interface = A_{\Gamma \Pi} tmp_vec_primal0
      MatMult(data->subSolver->getMatCoarse(1, 0), data->tmp_vec_primal0, data->tmp_vec_interface);
      // y_interface = A_{\Gamma B} tmp_vec_b0 + A_{\Gamma \Pi} tmp_vec_primal0
      VecAXPY(y_interface, 1.0, data->tmp_vec_interface);

      // y_lagrange = J tmp_vec_b0
      MatMult(*(data->mat_lagrange), data->tmp_vec_b0, y_lagrange);

      FetiTimings::fetiSolve += (MPI::Wtime() - wtime);


      return 0;
    }


    PetscErrorCode petscApplyFetiDirichletPrecon(PC pc, Vec x, Vec y)
    {
      double wtime = MPI::Wtime();

      // Get data for the preconditioner
      void* ctx;
      PCShellGetContext(pc, &ctx);
      FetiDirichletPreconData* data = static_cast<FetiDirichletPreconData*>(ctx);

      // Multiply with scaled Lagrange constraint matrix.
      MatMultTranspose(*(data->mat_lagrange_scaled), x, data->tmp_vec_b);


      // === Restriction of the B nodes to the boundary nodes. ===

      int nLocalB;
      int nLocalDuals;
      VecGetLocalSize(data->tmp_vec_b, &nLocalB);
      VecGetLocalSize(data->tmp_vec_duals0, &nLocalDuals);

      PetscScalar* local_b, *local_duals;
      VecGetArray(data->tmp_vec_b, &local_b);
      VecGetArray(data->tmp_vec_duals0, &local_duals);

      for (map<int, int>::iterator it = data->localToDualMap.begin();
           it != data->localToDualMap.end(); ++it)
        local_duals[it->second] = local_b[it->first];

      VecRestoreArray(data->tmp_vec_b, &local_b);
      VecRestoreArray(data->tmp_vec_duals0, &local_duals);


      // === K_DD - K_DI inv(K_II) K_ID ===

      MatMult(*(data->mat_duals_duals), data->tmp_vec_duals0, data->tmp_vec_duals1);

      MatMult(*(data->mat_interior_duals), data->tmp_vec_duals0, data->tmp_vec_interior);
      KSPSolve(*(data->ksp_interior), data->tmp_vec_interior, data->tmp_vec_interior);
      MatMult(*(data->mat_duals_interior), data->tmp_vec_interior, data->tmp_vec_duals0);

      VecAXPBY(data->tmp_vec_duals0, 1.0, -1.0, data->tmp_vec_duals1);


      // === Prolongation from local dual nodes to B nodes.

      VecGetArray(data->tmp_vec_b, &local_b);
      VecGetArray(data->tmp_vec_duals0, &local_duals);

      for (map<int, int>::iterator it = data->localToDualMap.begin();
           it != data->localToDualMap.end(); ++it)
        local_b[it->first] = local_duals[it->second];

      VecRestoreArray(data->tmp_vec_b, &local_b);
      VecRestoreArray(data->tmp_vec_duals0, &local_duals);


      // Multiply with scaled Lagrange constraint matrix.
      MatMult(*(data->mat_lagrange_scaled), data->tmp_vec_b, y);

      FetiTimings::fetiPreconditioner += (MPI::Wtime() - wtime);

      return 0;
    }


    PetscErrorCode petscApplyFetiLumpedPrecon(PC pc, Vec xvec, Vec yvec)
    {
      // Get data for the preconditioner
      void* ctx;
      PCShellGetContext(pc, &ctx);
      FetiLumpedPreconData* data = static_cast<FetiLumpedPreconData*>(ctx);

      // Multiply with scaled Lagrange constraint matrix.
      MatMultTranspose(*(data->mat_lagrange_scaled), xvec, data->tmp_vec_b0);


      // === Restriction of the B nodes to the boundary nodes. ===

      int nLocalB;
      int nLocalDuals;
      VecGetLocalSize(data->tmp_vec_b0, &nLocalB);
      VecGetLocalSize(data->tmp_vec_duals0, &nLocalDuals);

      PetscScalar* local_b, *local_duals;
      VecGetArray(data->tmp_vec_b0, &local_b);
      VecGetArray(data->tmp_vec_duals0, &local_duals);

      for (map<int, int>::iterator it = data->localToDualMap.begin();
           it != data->localToDualMap.end(); ++it)
        local_duals[it->second] = local_b[it->first];

      VecRestoreArray(data->tmp_vec_b0, &local_b);
      VecRestoreArray(data->tmp_vec_duals0, &local_duals);


      // === K_DD ===


      MatMult(*(data->mat_duals_duals), data->tmp_vec_duals0, data->tmp_vec_duals1);


      // === Prolongation from local dual nodes to B nodes.

      VecGetArray(data->tmp_vec_b0, &local_b);
      VecGetArray(data->tmp_vec_duals1, &local_duals);

      for (map<int, int>::iterator it = data->localToDualMap.begin();
           it != data->localToDualMap.end(); ++it)
        local_b[it->first] = local_duals[it->second];

      VecRestoreArray(data->tmp_vec_b0, &local_b);
      VecRestoreArray(data->tmp_vec_duals0, &local_duals);


      // Multiply with scaled Lagrange constraint matrix.
      MatMult(*(data->mat_lagrange_scaled), data->tmp_vec_b0, yvec);

      return 0;
    }


    PetscErrorCode petscApplyFetiInterfaceLumpedPrecon(PC pc, Vec xvec, Vec yvec)
    {
      FUNCNAME("precon");
      // Get data for the preconditioner
      void* ctx;
      PCShellGetContext(pc, &ctx);
      FetiInterfaceLumpedPreconData* data =
        static_cast<FetiInterfaceLumpedPreconData*>(ctx);

      Vec x_interface, x_lagrange, y_interface, y_lagrange;
      VecNestGetSubVec(xvec, 0, &x_interface);
      VecNestGetSubVec(xvec, 1, &x_lagrange);
      VecNestGetSubVec(yvec, 0, &y_interface);
      VecNestGetSubVec(yvec, 1, &y_lagrange);

      //VecCopy(x_interface, y_interface);
      KSPSolve(data->ksp_mass, x_interface, y_interface);

      MatMultTranspose(*(data->mat_lagrange_scaled), x_lagrange, data->tmp_vec_b0);

      // === Restriction of the B nodes to the boundary nodes. ===
      int nLocalB;
      int nLocalDuals;
      VecGetLocalSize(data->tmp_vec_b0, &nLocalB);
      VecGetLocalSize(data->tmp_vec_duals0, &nLocalDuals);

      PetscScalar* local_b, *local_duals;
      VecGetArray(data->tmp_vec_b0, &local_b);
      VecGetArray(data->tmp_vec_duals0, &local_duals);

      for (map<int, int>::iterator it = data->localToDualMap.begin();
           it != data->localToDualMap.end(); ++it)
        local_duals[it->second] = local_b[it->first];

      VecRestoreArray(data->tmp_vec_b0, &local_b);
      VecRestoreArray(data->tmp_vec_duals0, &local_duals);


      // === K_DD ===

      MatMult(*(data->mat_duals_duals), data->tmp_vec_duals0, data->tmp_vec_duals1);


      // === Prolongation from local dual nodes to B nodes.

      VecGetArray(data->tmp_vec_b0, &local_b);
      VecGetArray(data->tmp_vec_duals1, &local_duals);

      for (map<int, int>::iterator it = data->localToDualMap.begin();
           it != data->localToDualMap.end(); ++it)
        local_b[it->first] = local_duals[it->second];

      VecRestoreArray(data->tmp_vec_b0, &local_b);
      VecRestoreArray(data->tmp_vec_duals0, &local_duals);

      // Multiply with scaled Lagrange constraint matrix.
      MatMult(*(data->mat_lagrange_scaled), data->tmp_vec_b0, y_lagrange);

      return 0;
    }

  }
}
