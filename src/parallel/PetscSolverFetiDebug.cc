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

#include <mpi.h>
#include <petsc.h>
#include "parallel/ParallelDofMapping.h"
#include "parallel/PetscSolverFeti.h"
#include "parallel/PetscSolverFetiDebug.h"
#include "io/VtkWriter.h"

namespace AMDiS
{
  namespace Parallel
  {

    using namespace std;

    void PetscSolverFetiDebug::debugNullSpace(PetscSolverFeti& feti,
        SystemVector& vec)
    {
      FUNCNAME("PetscSolverFetiDebug::debugNullSpace()");

      TEST_EXIT(feti.stokesMode)
      ("This function works only for the stokes mode!\n");

      // === First, get the null space basis directly from the mat object and ===
      // === make a test whether this is a member of the null space.          ===

      {
        MatNullSpace nullSpace;
        MatGetNullSpace(feti.mat_feti, &nullSpace);

        PetscBool hasConst;
        PetscInt nBasisVec;
        const Vec* nullSpaceBasis;
        MatNullSpaceGetVecs(nullSpace, &hasConst, &nBasisVec, &nullSpaceBasis);

        TEST_EXIT(nBasisVec == 1)("Something is wrong!\n");

        Vec vecSol;
        VecDuplicate(nullSpaceBasis[0], &vecSol);
        MatMult(feti.mat_feti, nullSpaceBasis[0], vecSol);

        Vec vecSol0, vecSol1;
        VecNestGetSubVec(vecSol, 0, &vecSol0);
        VecNestGetSubVec(vecSol, 1, &vecSol1);

        PetscReal norm, norm0, norm1;
        VecNorm(vecSol, NORM_2, &norm);
        VecNorm(vecSol0, NORM_2, &norm0);
        VecNorm(vecSol1, NORM_2, &norm1);

        MSG("Null space norm: %e (%e %e)\n", norm, norm0, norm1);
        VecDestroy(&vecSol);
      }


      // === Try to reproduce the null space basis on the non reduced form of ===
      // === the FETI-DP system.                                              ===

      Mat fetiMat;
      createNestedFetiMat(feti, fetiMat);

      Vec vecArray[4];

      Vec ktest0, ktest1, ktest2;
      ParallelDofMapping& localDofMap = feti.localDofMap;
      feti.createLocalVec(localDofMap, ktest0);
      feti.createLocalVec(localDofMap, ktest1);
      feti.createVec(localDofMap, ktest2, feti.nGlobalOverallInterior);

      const FiniteElemSpace* pressureFeSpace =
        feti.componentSpaces[feti.pressureComponent];
      DofMap& m = localDofMap[feti.pressureComponent].getMap();
      for (DofMap::iterator it = m.begin(); it != m.end(); ++it)
      {
        if (feti.dofMap[pressureFeSpace].isRankDof(it->first))
        {
          int index = localDofMap.getLocalMatIndex(feti.pressureComponent, it->first);
          VecSetValue(ktest0, index, 1.0, INSERT_VALUES);
        }
      }
      VecAssemblyBegin(ktest0);
      VecAssemblyEnd(ktest0);
      MatMult(feti.subdomain->getMatInterior(), ktest0, ktest1);

      PetscScalar* valarray;
      VecGetArray(ktest0, &valarray);
      VecCreateMPIWithArray(PETSC_COMM_WORLD, 1,
                            localDofMap.getRankDofs(), feti.nGlobalOverallInterior,
                            valarray, &vecArray[0]);

      Vec ktest3;
      VecGetArray(ktest1, &valarray);
      VecCreateMPIWithArray(PETSC_COMM_WORLD, 1,
                            localDofMap.getRankDofs(), feti.nGlobalOverallInterior,
                            valarray, &ktest3);

      feti.createVec(feti.primalDofMap, vecArray[1]);
      VecSet(vecArray[1], 0.0);

      feti.createVec(feti.interfaceDofMap, vecArray[2]);
      VecSet(vecArray[2], 1.0);

      feti.createVec(feti.lagrangeMap, vecArray[3]);
      MatMult(feti.subdomain->getMatInteriorCoarse(1), vecArray[2], ktest2);
      VecAXPY(ktest2, 1.0, ktest3);
      MatMult(feti.mat_lagrange_scaled, ktest2, vecArray[3]);
      VecScale(vecArray[3], -1.0);


      Vec nullSpaceBasis;
      VecCreateNest(feti.domainComm, 4, PETSC_NULL, vecArray, &nullSpaceBasis);

      Vec vecSol;
      VecDuplicate(nullSpaceBasis, &vecSol);

      MatMult(fetiMat, nullSpaceBasis, vecSol);
      PetscReal norm;
      VecNorm(vecSol, NORM_2, &norm);
      MSG("Null space norm: %e\n", norm);

      Vec vec0;
      Vec vec1;
      Vec vec2;
      Vec vec3;
      VecNestGetSubVec(vecSol, 0, &vec0);
      VecNestGetSubVec(vecSol, 1, &vec1);
      VecNestGetSubVec(vecSol, 2, &vec2);
      VecNestGetSubVec(vecSol, 3, &vec3);

      PetscReal norm0, norm1, norm2, norm3;
      VecNorm(vec0, NORM_2, &norm0);
      VecNorm(vec1, NORM_2, &norm1);
      VecNorm(vec2, NORM_2, &norm2);
      VecNorm(vec3, NORM_2, &norm3);
      MSG("Splitted null space norm: %e %e %e %e\n", norm0, norm1, norm2, norm3);


      feti.recoverSolution(vec0, vec1, vec);
      feti.recoverInterfaceSolution(vec2, vec);

      io::VtkWriter::writeFile(&vec, "nullspace.vtu");

      MatDestroy(&fetiMat);
      VecDestroy(&nullSpaceBasis);
      VecDestroy(&vecSol);
    }


    void PetscSolverFetiDebug::createInteriorMat(PetscSolverFeti& feti, Mat& mat)
    {
      FUNCNAME("PetscSolverFetiDebug::createInteriorMat()");

      Mat& localMat = feti.subdomain->getMatInterior();
      int nRow, nCol;
      MatGetSize(localMat, &nRow, &nCol);

      TEST_EXIT_DBG(nRow == nCol)("Should not happen!\n");

      int* dnnz = new int[nRow];
      for (int i = 0; i < nRow; i++)
      {
        MatGetRow(localMat, i, &nCol, PETSC_NULL, PETSC_NULL);
        dnnz[i] = nCol;
        MatRestoreRow(localMat, i, &nCol, PETSC_NULL, PETSC_NULL);
      }

      MatCreateAIJ(PETSC_COMM_WORLD,
                   nRow, nRow,
                   feti.nGlobalOverallInterior, feti.nGlobalOverallInterior,
                   0, dnnz, 0, PETSC_NULL, &mat);

      PetscInt rStart, rEnd;
      MatGetOwnershipRange(mat, &rStart, &rEnd);

      for (int i = 0; i < nRow; i++)
      {
        const PetscInt* cols;
        const PetscScalar* vals;

        MatGetRow(localMat, i, &nCol, &cols, &vals);
        for (int j = 0; j < nCol; j++)
          MatSetValue(mat, i + rStart, cols[j] + rStart, vals[j], INSERT_VALUES);
      }

      MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

      delete [] dnnz;
    }


    void PetscSolverFetiDebug::createNestedFetiMat(PetscSolverFeti& feti, Mat& mat)
    {
      FUNCNAME("PetscSolverFetiDebug::createNestedFetiMat()");

      PetscSolver* subdomain = feti.subdomain;
      if (feti.stokesMode)
      {
        vector<Mat> nestMat;
        nestMat.resize(16);

        createInteriorMat(feti, nestMat[0]);
        nestMat[1] = subdomain->getMatInteriorCoarse(0);
        nestMat[2] = subdomain->getMatInteriorCoarse(1);
        MatTranspose(feti.mat_lagrange, MAT_INITIAL_MATRIX, &nestMat[3]);
        nestMat[4] = subdomain->getMatCoarseInterior(0);

        nestMat[5] = subdomain->getMatCoarse(0, 0);
        nestMat[6] = subdomain->getMatCoarse(0, 1);
        nestMat[7] = PETSC_NULL;
        nestMat[8] = subdomain->getMatCoarseInterior(1);
        nestMat[9] = subdomain->getMatCoarse(1, 0);
        nestMat[10] = PETSC_NULL;
        nestMat[11] = PETSC_NULL;
        nestMat[12] = feti.mat_lagrange;
        nestMat[13] = PETSC_NULL;
        nestMat[14] = PETSC_NULL;
        nestMat[15] = PETSC_NULL;

        //       Mat nestFetiMat;
        MatCreateNest(feti.domainComm, 4, PETSC_NULL, 4, PETSC_NULL,
                      &(nestMat[0]), &mat);
      }
      else
      {
        vector<Mat> nestMat;
        nestMat.resize(9);

        createInteriorMat(feti, nestMat[0]);
        nestMat[1] = subdomain->getMatInteriorCoarse(0);
        MatTranspose(feti.mat_lagrange, MAT_INITIAL_MATRIX, &nestMat[2]);
        nestMat[3] = subdomain->getMatCoarseInterior(0);
        nestMat[4] = subdomain->getMatCoarse(0, 0);
        nestMat[5] = PETSC_NULL;
        nestMat[6] = feti.mat_lagrange;
        nestMat[7] = PETSC_NULL;
        nestMat[8] = PETSC_NULL;

        //       Mat nestFetiMat;
        MatCreateNest(feti.domainComm, 3, PETSC_NULL, 3, PETSC_NULL,
                      &(nestMat[0]), &mat);
      }
    }


    void PetscSolverFetiDebug::createExplicitFetiMat(PetscSolverFeti& feti,
        Mat fetiMat,
        Mat& explicitMat,
        int& nnzCounter)
    {
      FUNCNAME("PetscSolverFetiDebug::createExplicitFetiMat()");

      ParallelDofMapping& lagrangeMap = feti.lagrangeMap;

      nnzCounter = 0;

      if (feti.stokesMode == false)
      {
        int nRankRows = lagrangeMap.getRankDofs();
        int nOverallRows = lagrangeMap.getOverallDofs();
        MatCreateAIJ(feti.domainComm,
                     nRankRows, nRankRows, nOverallRows, nOverallRows,
                     nOverallRows, PETSC_NULL, nOverallRows, PETSC_NULL,
                     &explicitMat);

        Vec unitVector;
        Vec resultVector;
        feti.createVec(lagrangeMap, unitVector);
        feti.createVec(lagrangeMap, resultVector);

        PetscInt low, high;
        VecGetOwnershipRange(unitVector, &low, &high);
        int nLocal = high - low;
        int nnzCounter = 0;

        for (int i = 0; i < lagrangeMap.getOverallDofs(); i++)
        {
          VecSet(unitVector, 0.0);
          if (i >= low && i < high)
            VecSetValue(unitVector, i, 1.0, INSERT_VALUES);
          VecAssemblyBegin(unitVector);
          VecAssemblyEnd(unitVector);

          MatMult(fetiMat, unitVector, resultVector);

          if (feti.fetiPreconditioner != FETI_NONE)
          {
            PCApply(feti.precon_feti, resultVector, unitVector);
            VecCopy(unitVector, resultVector);
          }

          PetscScalar* vals;
          VecGetArray(resultVector, &vals);
          for (int j = 0; j < nLocal; j++)
          {
            if (fabs(vals[j]) > 1e-30)
            {
              MatSetValue(explicitMat, low + j, i, vals[j], INSERT_VALUES);
              nnzCounter++;
            }
          }
          VecRestoreArray(resultVector, &vals);
        }

        VecDestroy(&unitVector);
        VecDestroy(&resultVector);
      }
      else
      {
        ParallelDofMapping& interfaceDofMap = feti.interfaceDofMap;

        int nRankRows = interfaceDofMap.getRankDofs() + lagrangeMap.getRankDofs();
        int nOverallRows =
          interfaceDofMap.getOverallDofs() + lagrangeMap.getOverallDofs();

        MatCreateAIJ(feti.domainComm,
                     nRankRows, nRankRows, nOverallRows, nOverallRows,
                     nOverallRows, PETSC_NULL, nOverallRows, PETSC_NULL,
                     &explicitMat);

        Vec unitVector[2];
        Vec resultVector[2];
        feti.createVec(interfaceDofMap, unitVector[0]);
        feti.createVec(interfaceDofMap, resultVector[0]);
        feti.createVec(lagrangeMap, unitVector[1]);
        feti.createVec(lagrangeMap, resultVector[1]);

        Vec unitNestVec, resultNestVec;
        VecCreateNest(feti.domainComm, 2, PETSC_NULL, unitVector, &unitNestVec);
        VecCreateNest(feti.domainComm, 2, PETSC_NULL, resultVector, &resultNestVec);

        PetscInt lowInterface, highInterface;
        VecGetOwnershipRange(unitVector[0], &lowInterface, &highInterface);
        int nLocalInterface = highInterface - lowInterface;

        PetscInt lowLagrange, highLagrange;
        VecGetOwnershipRange(unitVector[1], &lowLagrange, &highLagrange);
        int nLocalLagrange = highLagrange - lowLagrange;

        VecSet(unitVector[1], 0.0);
        for (int i = 0; i < interfaceDofMap.getOverallDofs(); i++)
        {
          VecSet(unitVector[0], 0.0);
          if (i >= lowInterface && i < highInterface)
            VecSetValue(unitVector[0], i, 1.0, INSERT_VALUES);
          VecAssemblyBegin(unitVector[0]);
          VecAssemblyEnd(unitVector[0]);

          VecAssemblyBegin(unitNestVec);
          VecAssemblyEnd(unitNestVec);

          MatMult(fetiMat, unitNestVec, resultNestVec);

          PetscScalar* vals;
          VecGetArray(resultVector[0], &vals);
          for (int j = 0; j < nLocalInterface; j++)
          {
            if (fabs(vals[j]) > 1e-30)
            {
              MatSetValue(explicitMat, lowInterface + j, i, vals[j], INSERT_VALUES);
              nnzCounter++;
            }
          }
          VecRestoreArray(resultVector[0], &vals);

          VecGetArray(resultVector[1], &vals);
          for (int j = 0; j < nLocalLagrange; j++)
          {
            if (fabs(vals[j]) > 1e-30)
            {
              MatSetValue(explicitMat,
                          interfaceDofMap.getOverallDofs() + lowLagrange + j, i,
                          vals[j], INSERT_VALUES);
              nnzCounter++;
            }
          }
          VecRestoreArray(resultVector[1], &vals);
        }


        VecSet(unitVector[0], 0.0);
        for (int i = 0; i < lagrangeMap.getOverallDofs(); i++)
        {
          VecSet(unitVector[1], 0.0);
          if (i >= lowLagrange && i < highLagrange)
            VecSetValue(unitVector[1], i, 1.0, INSERT_VALUES);
          VecAssemblyBegin(unitVector[1]);
          VecAssemblyEnd(unitVector[1]);

          VecAssemblyBegin(unitNestVec);
          VecAssemblyEnd(unitNestVec);


          MatMult(fetiMat, unitNestVec, resultNestVec);

          PetscScalar* vals;
          VecGetArray(resultVector[0], &vals);
          for (int j = 0; j < nLocalInterface; j++)
          {
            if (fabs(vals[j]) > 1e-30)
            {
              MatSetValue(explicitMat,
                          lowInterface + j,
                          interfaceDofMap.getOverallDofs() + i,
                          vals[j], INSERT_VALUES);
              nnzCounter++;
            }
          }
          VecRestoreArray(resultVector[0], &vals);

          VecGetArray(resultVector[1], &vals);
          for (int j = 0; j < nLocalLagrange; j++)
          {
            if (fabs(vals[j]) > 1e-30)
            {
              MatSetValue(explicitMat,
                          interfaceDofMap.getOverallDofs() + lowLagrange + j,
                          interfaceDofMap.getOverallDofs() + i,
                          vals[j], INSERT_VALUES);
              nnzCounter++;
            }
          }
          VecRestoreArray(resultVector[1], &vals);
        }

        VecDestroy(&unitVector[0]);
        VecDestroy(&unitVector[1]);
        VecDestroy(&resultVector[0]);
        VecDestroy(&resultVector[1]);
        VecDestroy(&unitNestVec);
        VecDestroy(&resultNestVec);
      }

      MatAssemblyBegin(explicitMat, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(explicitMat, MAT_FINAL_ASSEMBLY);
      mpi::globalAdd(nnzCounter);
    }


    void PetscSolverFetiDebug::createExplicitVec(PetscSolverFeti& feti,
        Vec nestedVec,
        Vec& explicitVec)
    {
      FUNCNAME("PetscSolverFetiDebug::createExplicitVec()");

      int nNested = 0;
      VecNestGetSize(nestedVec, &nNested);

      TEST_EXIT_DBG(nNested == 2)
      ("Only supported for 2-nest vectors, not for %d-nest vectors!\n", nNested);

      Vec v0, v1;
      VecNestGetSubVec(nestedVec, 0, &v0);
      VecNestGetSubVec(nestedVec, 1, &v1);

      int s0, s1;
      VecGetSize(v0, &s0);
      VecGetSize(v1, &s1);

      int l0, l1;
      VecGetLocalSize(v0, &l0);
      VecGetLocalSize(v1, &l1);

      int rStart0, rStart1;
      VecGetOwnershipRange(v0, &rStart0, PETSC_NULL);
      VecGetOwnershipRange(v1, &rStart1, PETSC_NULL);

      VecCreateMPI(feti.domainComm, l0 + l1, s0 + s1, &explicitVec);

      double* val;
      VecGetArray(v0, &val);
      for (int i = 0; i < l0; i++)
        VecSetValue(explicitVec, rStart0 + i, val[i], INSERT_VALUES);
      VecRestoreArray(v0, &val);

      VecGetArray(v1, &val);
      for (int i = 0; i < l1; i++)
        VecSetValue(explicitVec, s0 + rStart1 + i, val[i], INSERT_VALUES);
      VecRestoreArray(v1, &val);

      VecAssemblyBegin(explicitVec);
      VecAssemblyEnd(explicitVec);
    }


    void PetscSolverFetiDebug::writeNullSpace(PetscSolverFeti& feti,
        Vec nullSpaceBasis)
    {
      FUNCNAME("PetscSolverFetiDebug::writeNullSpace()");

      int writeFetiSystem = 0;
      Parameters::get("parallel->debug->write feti system",
                      writeFetiSystem);
      if (writeFetiSystem)
      {
        Vec vecTmp;
        createExplicitVec(feti, nullSpaceBasis, vecTmp);
        PetscViewer petscView;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "nullspace.vec",
                              FILE_MODE_WRITE, &petscView);
        VecView(vecTmp, petscView);
        PetscViewerDestroy(&petscView);

        VecDestroy(&vecTmp);
        MSG("Written FETI-DP null space basis vector!\n");
      }
    }


    void PetscSolverFetiDebug::debugFeti(PetscSolverFeti& feti, Vec dbgRhsVec)
    {
      FUNCNAME("PetscSolverFetiDebug:::debugFeti()");

      int writeFetiSystem = 0;
      Parameters::get("parallel->debug->write feti system",
                      writeFetiSystem);
      if (!writeFetiSystem)
        return;

      MSG("Start creating explicit FETI-DP matrix!\n");

      Mat fetiMat;
      int nnzCounter = 0;
      createExplicitFetiMat(feti, feti.mat_feti, fetiMat, nnzCounter);

      PetscViewer petscView;
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, "feti.mat",
                            FILE_MODE_WRITE, &petscView);
      MatView(fetiMat, petscView);
      PetscViewerDestroy(&petscView);

      bool a = feti.testMatrixSymmetric(fetiMat, true);
      MSG("SYMMETRIC TEST: %d\n", a);

      int r, c;
      MatGetSize(fetiMat, &r, &c);

      MSG("FETI-DP matrix written: %d x %d mat with %d nnz\n",
          r, c, nnzCounter);

      MatDestroy(&fetiMat);

      Vec rhsVec;
      createExplicitVec(feti, dbgRhsVec, rhsVec);

      PetscViewerBinaryOpen(PETSC_COMM_WORLD, "feti_rhs.vec",
                            FILE_MODE_WRITE, &petscView);
      VecView(rhsVec, petscView);
      PetscViewerDestroy(&petscView);

      VecDestroy(&rhsVec);
    }


    int PetscSolverFetiDebug::testZeroRows(Mat mat)
    {
      FUNCNAME("PetscSolverFetiDebug::testZeroRows()");

      int nZeroRows = 0;

      int firstIndex, lastIndex;
      MatGetOwnershipRange(mat, &firstIndex, &lastIndex);

      for (int i = firstIndex; i < lastIndex; i++)
      {
        PetscInt nCols;
        const PetscScalar* vals;

        MatGetRow(mat, i, &nCols, PETSC_NULL, &vals);
        bool nonZeroValues = false;
        for (int j = 0; j < nCols; j++)
        {
          if (vals[j] != 0.0)
          {
            nonZeroValues = true;
            break;
          }
        }

        if (!nonZeroValues)
        {
          MSG("Zero matrix row for global row index %d\n", i);
          nZeroRows++;
        }

        MatRestoreRow(mat, i, &nCols, PETSC_NULL, &vals);
      }

      return nZeroRows;
    }


    void PetscSolverFetiDebug::writePrimalFiles(PetscSolverFeti& feti)
    {
      FUNCNAME("PetscSolverFetiDebug::writePrimalFiles()");

      Mesh* mesh = feti.componentSpaces[0]->getMesh();
      DOFVector<WorldVector<double>> coords(feti.componentSpaces[0], "coords");
      mesh->getDofIndexCoords(coords);


      // === Write vertex primals. ===

      {
        StdMpi<vector<double>> stdMpi(MPI::COMM_WORLD);
        vector<double> sendData;

        DofMap& primals = feti.primalDofMap[0].getMap();
        for (DofMap::iterator it = primals.begin(); it != primals.end(); ++it)
        {
          if (feti.primalDofMap[0].isRankDof(it->first))
          {
            for (int i = 0; i < mesh->getDim(); i++)
              sendData.push_back(coords[it->first][i]);
          }
        }

        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
          for (int i = 1; i < MPI::COMM_WORLD.Get_size(); i++)
            stdMpi.recv(i);
        }
        else
        {
          stdMpi.send(0, sendData);
        }

        stdMpi.startCommunication();

        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
          vector<WorldVector<double>> allPrimals;

          for (map<int, vector<double>>::iterator it = stdMpi.getRecvData().begin();
               it != stdMpi.getRecvData().end(); ++it)
          {
            vector<double>& recvData = it->second;

            TEST_EXIT(recvData.size() % mesh->getDim() == 0)
            ("Wrong number of coordinates received!\n");

            int nCoords = recvData.size() / mesh->getDim();
            int counter = 0;
            for (int j = 0; j < nCoords; j++)
            {
              WorldVector<double> c;
              for (int k = 0; k < mesh->getDim(); k++)
                c[k] = recvData[counter++];
              allPrimals.push_back(c);
            }
          }

          ofstream file;
          file.open("primals.xyz");
          file << allPrimals.size() << "\n";
          file << "Primals\n";
          for (int i = 0; i < static_cast<int>(allPrimals.size()); i++)
          {
            file << "P ";
            if (mesh->getDim() == 2)
              file << allPrimals[i][0] << " " << allPrimals[i][1] << " 0.0";
            else
              file << allPrimals[i][0] << " " << allPrimals[i][1] << " " << allPrimals[i][2];
            file << "\n";
          }
          file.close();
        }
      }


      // === Write face primals. ===

      if (mesh->getDim() == 3)
      {
        StdMpi<vector<double>> stdMpi(MPI::COMM_WORLD);
        vector<double> sendData;

        vector<vector<BoundaryObject>> coarseSpace = feti.getCoarseFaces();
        for (vector<vector<BoundaryObject>>::iterator it = coarseSpace.begin();
             it != coarseSpace.end(); ++it)
        {
          for (vector<BoundaryObject>::iterator edgeIt = it->begin();
               edgeIt != it->end(); ++edgeIt)
          {
            if (edgeIt->subObj == FACE)
            {
              DofFace face = edgeIt->el->getFace(edgeIt->ithObj);
              WorldVector<double> c0 = coords[std::get<0>(face)];
              WorldVector<double> c1 = coords[std::get<1>(face)];
              WorldVector<double> c2 = coords[std::get<2>(face)];
              sendData.push_back(c0[0]);
              sendData.push_back(c0[1]);
              sendData.push_back(c0[2]);
              sendData.push_back(c1[0]);
              sendData.push_back(c1[1]);
              sendData.push_back(c1[2]);
              sendData.push_back(c2[0]);
              sendData.push_back(c2[1]);
              sendData.push_back(c2[2]);
            }
          }
        }

        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
          TEST_EXIT(sendData.size() == 0)("Should not happen!\n");
        }

        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
          for (int i = 1; i < MPI::COMM_WORLD.Get_size(); i++)
            stdMpi.recv(i);
        }
        else
        {
          stdMpi.send(0, sendData);
        }

        stdMpi.startCommunication();

        if (MPI::COMM_WORLD.Get_rank() == 0)
        {
          vector<WorldVector<double>> faceNodes;

          for (map<int, vector<double>>::iterator it = stdMpi.getRecvData().begin();
               it != stdMpi.getRecvData().end(); ++it)
          {
            vector<double>& recvData = it->second;

            TEST_EXIT(recvData.size() % 9 == 0)
            ("Wrong number of coordinates received!\n");

            int nNodes = recvData.size() / 3;
            int counter = 0;
            for (int j = 0; j < nNodes; j++)
            {
              WorldVector<double> c;
              c[0] = recvData[counter++];
              c[1] = recvData[counter++];
              c[2] = recvData[counter++];
              faceNodes.push_back(c);
            }
          }

          TEST_EXIT(faceNodes.size() % 3 == 0)("Should not happen!\n");

          int nElements = faceNodes.size() / 3;

          ofstream file;
          file.open("faces.vtu");

          file << "<?xml version=\"1.0\"?>\n";
          file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
          file << "  <UnstructuredGrid>\n";
          file << "    <Piece NumberOfPoints=\"" << faceNodes.size()
               << "\" NumberOfCells=\"" << nElements << "\">\n";
          file << "      <Points>\n";
          file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

          for (int i = 0; i < static_cast<int>(faceNodes.size()); i++)
            file << " " << faceNodes[i][0] << " " << faceNodes[i][1] << " " << faceNodes[i][2] << "\n";

          file << "        </DataArray>\n";
          file << "      </Points>\n";
          file << "      <Cells>\n";
          file << "        <DataArray type=\"Int32\" Name=\"offsets\">\n";

          for (int i = 0; i < nElements; i++)
            file << " " << (i + 1) * 3 << "\n";

          file << "        </DataArray>\n";
          file << "        <DataArray type=\"UInt8\" Name=\"types\">\n";

          for (int i = 0; i < nElements; i++)
            file << " 5\n";

          file << "        </DataArray>\n";
          file << "        <DataArray type=\"Int32\" Name=\"connectivity\">\n";

          {
            int counter = 0;
            for (int i = 0; i < nElements; i++)
            {
              file << " " << counter << " " << counter + 1 << " " << counter + 2 << "\n";
              counter += 3;
            }
          }

          file << "        </DataArray>\n";
          file << "      </Cells>\n";
          file << "      <PointData>\n";

          file << "        <DataArray type=\"Float32\" Name=\"value\" format=\"ascii\">\n";

          for (int i = 0; i < static_cast<int>(faceNodes.size()); i++)
            file << " 0.0\n";

          file << "        </DataArray>\n";

          file << "      </PointData>\n";
          file << "    </Piece>\n";
          file << "  </UnstructuredGrid>\n";
          file << "</VTKFile>\n";

          file.close();
        }
      }

    }

  }
}
