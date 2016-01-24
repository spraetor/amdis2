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


#include "parallel/PetscSolverSchur.h"
#include "parallel/StdMpi.h"
#include "parallel/MpiHelper.h"

namespace AMDiS
{
  namespace Parallel
  {

    using namespace std;

    void PetscSolverSchur::updateDofData(int nComponents)
    {
      FUNCNAME("PetscSolverSchur::updateDofData()");

      TEST_EXIT_DBG(meshDistributor)("No mesh distributor object defined!\n");
      TEST_EXIT_DBG(interiorMap)("No parallel DOF map defined!\n");

      const FiniteElemSpace* feSpace = componentSpaces[0];
      Mesh* mesh = feSpace->getMesh();
      //     typedef map<int, DofContainer> RankToDofContainer;
#if 0
      typedef map<DegreeOfFreedom, bool> DofIndexToBool;
#endif

      boundaryDofs.clear();
      std::set<DegreeOfFreedom> boundaryLocalDofs;
      for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, 0).getSendDofs(), feSpace);
           !it.end(); it.nextRank())
        for (; !it.endDofIter(); it.nextDof())
        {
          boundaryLocalDofs.insert(it.getDofIndex());
          boundaryDofs.insert((*interiorMap)[feSpace][it.getDofIndex()].global);
        }


      nBoundaryDofs = boundaryDofs.size();
      mpi::getDofNumbering(domainComm, nBoundaryDofs,
                           rStartBoundaryDofs, nOverallBoundaryDofs);


      DofContainerSet& edgeDofs =
        meshDistributor->getBoundaryDofInfo(feSpace, 0).geoDofs[EDGE];
      DofContainerSet& vertexDofs =
        meshDistributor->getBoundaryDofInfo(feSpace, 0).geoDofs[VERTEX];
      int nEdgeDofs = edgeDofs.size();
      int nVertexDofs = vertexDofs.size();

      TEST_EXIT_DBG(nEdgeDofs + nVertexDofs == nBoundaryDofs)
      ("Should not happen!\n");

      int rStartEdgeDofs, nOverallEdgeDofs;
      mpi::getDofNumbering(domainComm, nEdgeDofs,
                           rStartEdgeDofs, nOverallEdgeDofs);

      int rStartVertexDofs, nOverallVertexDofs;
      mpi::getDofNumbering(domainComm, nVertexDofs,
                           rStartVertexDofs, nOverallVertexDofs);

      TEST_EXIT_DBG(nOverallEdgeDofs + nOverallVertexDofs == nOverallBoundaryDofs)
      ("Should not happen!\n");


      mapGlobalBoundaryDof.clear();
#if 1
      {
        int counter = rStartEdgeDofs;
        for (DofContainerSet::iterator it = edgeDofs.begin();
             it != edgeDofs.end(); ++it)
          mapGlobalBoundaryDof[(*interiorMap)[feSpace][** it].global] =
            counter++;
      }
      {
        int counter = nOverallEdgeDofs + rStartVertexDofs;
        for (DofContainerSet::iterator it = vertexDofs.begin();
             it != vertexDofs.end(); ++it)
          mapGlobalBoundaryDof[(*interiorMap)[feSpace][** it].global] =
            counter++;
      }
#else
      {
        int counter = rStartBoundaryDofs;
        for (std::set<DegreeOfFreedom>::iterator it = boundaryDofs.begin();
             it != boundaryDofs.end(); ++it)
          mapGlobalBoundaryDof[*it] = counter++;
      }
#endif


      std::set<DegreeOfFreedom> otherBoundaryLocalDofs;
      for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, 0).getRecvDofs(), feSpace);
           !it.end(); it.nextRank())
        for (; !it.endDofIter(); it.nextDof())
          otherBoundaryLocalDofs.insert(it.getDofIndex());

      interiorDofs.clear();

      ERROR_EXIT("Rewrite the following code block!\n");
#if 0
      DofIndexToBool& isRankDof = meshDistributor->getIsRankDof(feSpace);
      for (DofIndexToBool::iterator dofIt = isRankDof.begin();
           dofIt != isRankDof.end(); ++dofIt)
      {
        if (dofIt->second &&
            boundaryLocalDofs.count(dofIt->first) == 0 &&
            otherBoundaryLocalDofs.count(dofIt->first) == 0)
          interiorDofs.insert(meshDistributor->mapDofToGlobal(feSpace, dofIt->first));
      }
#endif

      nInteriorDofs = interiorDofs.size();
      mpi::getDofNumbering(domainComm, nInteriorDofs,
                           rStartInteriorDofs, nOverallInteriorDofs);

      {
        int counter = rStartInteriorDofs;
        mapGlobalInteriorDof.clear();
        for (std::set<DegreeOfFreedom>::iterator it = interiorDofs.begin();
             it != interiorDofs.end(); ++it)
          mapGlobalInteriorDof[*it] = counter++;
      }


      TEST_EXIT_DBG(nInteriorDofs > 0)("Should not happen!\n");


      StdMpi<vector<DegreeOfFreedom>> stdMpi(domainComm);
      for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, 0).getSendDofs(), feSpace);
           !it.end(); it.nextRank())
      {
        stdMpi.getSendData(it.getRank()).resize(0);
        stdMpi.getSendData(it.getRank()).reserve(it.getDofs().size());

        for (; !it.endDofIter(); it.nextDof())
        {
          int globalSendDof = (*interiorMap)[feSpace][it.getDofIndex()].global;

          TEST_EXIT_DBG(mapGlobalBoundaryDof.count(globalSendDof))
          ("No mapping for boundary DOF %d!\n", globalSendDof);

          stdMpi.getSendData(it.getRank()).push_back(mapGlobalBoundaryDof[globalSendDof]);
        }
      }

      stdMpi.updateSendDataSize();

      for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, 0).getRecvDofs(), feSpace);
           !it.end(); it.nextRank())
        stdMpi.recv(it.getRank());

      stdMpi.startCommunication();

      for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, 0).getRecvDofs(), feSpace);
           !it.end(); it.nextRank())
        for (; !it.endDofIter(); it.nextDof())
        {
          int globalRecvDof = (*interiorMap)[feSpace][it.getDofIndex()].global;
          mapGlobalBoundaryDof[globalRecvDof] =
            stdMpi.getRecvData(it.getRank())[it.getDofCounter()];
          boundaryDofs.insert(globalRecvDof);
        }


      // === Create PETSc IS structurs for interior and boundary DOFs. ===

      ISCreateStride(domainComm,
                     nInteriorDofs * nComponents,
                     (rStartInteriorDofs + rStartBoundaryDofs) * nComponents,
                     1, &interiorIs);

      ISCreateStride(domainComm,
                     nBoundaryDofs * nComponents,
                     (rStartInteriorDofs + rStartBoundaryDofs + nInteriorDofs) * nComponents,
                     1, &boundaryIs);
    }


    void PetscSolverSchur::fillPetscMatrix(Matrix<DOFMatrix*>* seqMat)
    {
      //     FUNCNAME("PetscSolverSchur::fillPetscMatrix()");

      createMatVec(*seqMat);

      const FiniteElemSpace* feSpace = componentSpaces[0];
      int nComponents = seqMat->getNumRows();
      updateDofData(nComponents);

      int nInteriorRows = nInteriorDofs * nComponents;
      int nOverallInteriorRows = nOverallInteriorDofs * nComponents;
      int nBoundaryRows = nBoundaryDofs * nComponents;
      int nOverallBoundaryRows = nOverallBoundaryDofs * nComponents;


      MatCreateAIJ(domainComm,
                   nInteriorRows, nInteriorRows,
                   nOverallInteriorRows, nOverallInteriorRows,
                   100, PETSC_NULL, 100, PETSC_NULL, &matA11);

      MatCreateAIJ(domainComm,
                   nBoundaryRows, nBoundaryRows,
                   nOverallBoundaryRows, nOverallBoundaryRows,
                   100, PETSC_NULL, 100, PETSC_NULL, &matA22);

      MatCreateAIJ(domainComm,
                   nInteriorRows, nBoundaryRows,
                   nOverallInteriorRows, nOverallBoundaryRows,
                   100, PETSC_NULL, 100, PETSC_NULL, &matA12);

      MatCreateAIJ(domainComm,
                   nBoundaryRows, nInteriorRows,
                   nOverallBoundaryRows, nOverallInteriorRows,
                   100, PETSC_NULL, 100, PETSC_NULL, &matA21);


      for (int i = 0; i < nComponents; i++)
        for (int j = 0; j < nComponents; j++)
          if ((*seqMat)[i][j])
            setDofMatrix((*seqMat)[i][j], nComponents, i, j);

      MatAssemblyBegin(matA11, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(matA11, MAT_FINAL_ASSEMBLY);

      MatAssemblyBegin(matA12, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(matA12, MAT_FINAL_ASSEMBLY);

      MatAssemblyBegin(matA21, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(matA21, MAT_FINAL_ASSEMBLY);

      MatAssemblyBegin(matA22, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(matA22, MAT_FINAL_ASSEMBLY);

      Mat tmpMat[2][2];
      tmpMat[0][0] = matA11;
      tmpMat[0][1] = matA12;
      tmpMat[1][0] = matA21;
      tmpMat[1][1] = matA22;

      IS tmpIS[2];
      tmpIS[0] = interiorIs;
      tmpIS[1] = boundaryIs;

      MatCreateNest(domainComm, 2, &tmpIS[0], 2, &tmpIS[0], &tmpMat[0][0],
                    &getMatInterior());
      MatNestSetVecType(getMatInterior(), VECNEST);

      matAssembly();

      int nRankRows = (*interiorMap)[feSpace].nRankDofs * nComponents;
      int nOverallRows = (*interiorMap)[feSpace].nOverallDofs * nComponents;

      VecCreateMPI(domainComm, nRankRows, nOverallRows, &petscSolVec);
    }


    void PetscSolverSchur::fillPetscRhs(SystemVector* vec)
    {
      //     FUNCNAME("PetscSolverSchur::fillPetscRhs()");

      const FiniteElemSpace* feSpace = componentSpaces[0];
      int nComponents = vec->getSize();
      int nRankRows = (*interiorMap)[feSpace].nRankDofs * nComponents;
      int nOverallRows = (*interiorMap)[feSpace].nOverallDofs * nComponents;

      VecCreateMPI(domainComm, nRankRows, nOverallRows, &(getVecRhsInterior()));

      for (int i = 0; i < nComponents; i++)
        setDofVector(getVecRhsInterior(), vec->getDOFVector(i), nComponents, i);

      vecRhsAssembly();
    }


    void PetscSolverSchur::solvePetscMatrix(SystemVector& vec,
                                            AdaptInfo& adaptInfo)
    {
      //     FUNCNAME("PetscSolverSchur::solvePetscMatrix()");

      const FiniteElemSpace* feSpace = componentSpaces[0];
      int nComponents = vec.getSize();

      KSPCreate(domainComm, &kspInterior);

#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(kspInterior, getMatInterior(), getMatInterior());
#else
      KSPSetOperators(kspInterior, getMatInterior(), getMatInterior(), SAME_NONZERO_PATTERN);
#endif
      KSPSetTolerances(kspInterior, 0.0, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT);
      KSPSetFromOptions(kspInterior);

      KSPGetPC(kspInterior, &pcInterior);
      PCSetType(pcInterior, PCFIELDSPLIT);
      PCFieldSplitSetIS(pcInterior, "interior", interiorIs);
      PCFieldSplitSetIS(pcInterior, "boundary", boundaryIs);
      PCSetFromOptions(pcInterior);


      KSPSolve(kspInterior, getVecRhsInterior(), petscSolVec);

      // === Transfere values from PETSc's solution vectors to AMDiS vectors. ===

      PetscScalar* vecPointer;
      VecGetArray(petscSolVec, &vecPointer);

      for (int i = 0; i < nComponents; i++)
      {
        DOFVector<double>::Iterator dofIt(vec.getDOFVector(i), USED_DOFS);
        for (dofIt.reset(); !dofIt.end(); ++dofIt)
        {
          DegreeOfFreedom globalRowDof = (*interiorMap)[feSpace][dofIt.getDOFIndex()].global;
          if (boundaryDofs.count(globalRowDof))
          {
            int index =
              (mapGlobalBoundaryDof[globalRowDof] - rStartBoundaryDofs + nInteriorDofs) * (i + 1);
            *dofIt = vecPointer[index];
          }
          else
          {
            int index =
              (mapGlobalInteriorDof[globalRowDof] - rStartInteriorDofs) * (i + 1);
            *dofIt = vecPointer[index];
          }
        }
      }

      VecRestoreArray(petscSolVec, &vecPointer);


      // === Synchronize DOFs at common DOFs, i.e., DOFs that correspond to ===
      // === more than one partition.                                       ===
      meshDistributor->synchVector(vec);


      // === Destroy PETSC's variables. ===

      MatDestroy(&matA11);
      MatDestroy(&matA12);
      MatDestroy(&matA21);
      MatDestroy(&matA22);

      matDestroy();
      vecDestroy();

      KSPDestroy(&kspInterior);
    }


    void PetscSolverSchur::setDofMatrix(DOFMatrix* seqMat, int dispMult,
                                        int dispAddRow, int dispAddCol)
    {
      FUNCNAME("PetscSolverSchur::setDofMatrix()");

      TEST_EXIT(seqMat)("No DOFMatrix!\n");

      const FiniteElemSpace* feSpace = componentSpaces[0];
      using mtl::tag::row;
      using mtl::tag::nz;
      using mtl::begin;
      using mtl::end;
      namespace traits= mtl::traits;
      typedef DOFMatrix::base_matrix_type Matrix;

      traits::col<Matrix>::type col(seqMat->getBaseMatrix());
      traits::const_value<Matrix>::type value(seqMat->getBaseMatrix());

      typedef traits::range_generator<row, Matrix>::type cursor_type;
      typedef traits::range_generator<nz, cursor_type>::type icursor_type;

      vector<int> colsBoundary, colsInterior;
      vector<double> valuesBoundary, valuesInterior;
      colsBoundary.reserve(300);
      colsInterior.reserve(300);
      valuesBoundary.reserve(300);
      valuesInterior.reserve(300);

      for (cursor_type cursor = begin<row>(seqMat->getBaseMatrix()),
           cend = end<row>(seqMat->getBaseMatrix()); cursor != cend; ++cursor)
      {

        // Global index of the current row DOF.
        int globalRowDof = (*interiorMap)[feSpace][cursor.value()].global;

        colsBoundary.clear();
        colsInterior.clear();
        valuesBoundary.clear();
        valuesInterior.clear();

        for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor);
             icursor != icend; ++icursor)
        {
          int globalColDof = (*interiorMap)[feSpace][col(*icursor)].global;

          if (boundaryDofs.count(globalColDof))
          {
            TEST_EXIT_DBG(mapGlobalBoundaryDof.count(globalColDof))
            ("Should not happen!\n");

            int colIndex =
              mapGlobalBoundaryDof[globalColDof] * dispMult + dispAddCol;

            colsBoundary.push_back(colIndex);
            valuesBoundary.push_back(value(*icursor));
          }
          else
          {
            TEST_EXIT_DBG(mapGlobalInteriorDof.count(globalColDof))
            ("Cannot find global interior mapping for global column DOF %d!\n",
             globalColDof);

            int colIndex =
              mapGlobalInteriorDof[globalColDof] * dispMult + dispAddCol;

            colsInterior.push_back(colIndex);
            valuesInterior.push_back(value(*icursor));
          }
        }

        if (boundaryDofs.count(globalRowDof))
        {
          TEST_EXIT_DBG(mapGlobalBoundaryDof.count(globalRowDof))
          ("Should not happen!\n");

          int rowIndex =
            mapGlobalBoundaryDof[globalRowDof] * dispMult + dispAddRow;

          MatSetValues(matA22, 1, &rowIndex, colsBoundary.size(),
                       &(colsBoundary[0]), &(valuesBoundary[0]), ADD_VALUES);
          MatSetValues(matA21, 1, &rowIndex, colsInterior.size(),
                       &(colsInterior[0]), &(valuesInterior[0]), ADD_VALUES);
        }
        else
        {
          TEST_EXIT_DBG(mapGlobalInteriorDof.count(globalRowDof))
          ("Cannot find global interior mapping for global row DOF %d!\n",
           globalRowDof);

          int rowIndex =
            mapGlobalInteriorDof[globalRowDof] * dispMult + dispAddRow;

          MatSetValues(matA11, 1, &rowIndex, colsInterior.size(),
                       &(colsInterior[0]), &(valuesInterior[0]), ADD_VALUES);
          MatSetValues(matA12, 1, &rowIndex, colsBoundary.size(),
                       &(colsBoundary[0]), &(valuesBoundary[0]), ADD_VALUES);
        }
      }
    }


    void PetscSolverSchur::setDofVector(Vec& petscVec, DOFVector<double>* vec,
                                        int dispMult, int dispAdd, bool rankOnly)
    {
      FUNCNAME_DBG("PetscSolverSchur::setDofVector()");

      const FiniteElemSpace* feSpace = componentSpaces[0];

      DOFVector<double>::Iterator dofIt(vec, USED_DOFS);
      for (dofIt.reset(); !dofIt.end(); ++dofIt)
      {
        if (rankOnly && !(*interiorMap)[feSpace].isRankDof(dofIt.getDOFIndex()))
          continue;

        // Calculate global row index of the DOF.
        DegreeOfFreedom globalRowDof =
          (*interiorMap)[feSpace][dofIt.getDOFIndex()].global;
        double value = *dofIt;

        if (boundaryDofs.count(globalRowDof))
        {
          TEST_EXIT_DBG(mapGlobalBoundaryDof.count(globalRowDof))
          ("Should not happen!\n");

          int index =
            (rStartInteriorDofs +
             nInteriorDofs +
             mapGlobalBoundaryDof[globalRowDof]) * dispMult + dispAdd;
          VecSetValues(petscVec, 1, &index, &value, ADD_VALUES);
        }
        else
        {
          TEST_EXIT_DBG(mapGlobalInteriorDof.count(globalRowDof))
          ("Should not happen!\n");

          int index =
            (rStartBoundaryDofs +
             mapGlobalInteriorDof[globalRowDof]) * dispMult + dispAdd;
          VecSetValues(petscVec, 1, &index, &value, ADD_VALUES);
        }
      }
    }

  }
} // end namespace Parallel, AMDiS
