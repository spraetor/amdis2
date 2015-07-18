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


#include "parallel/ParallelCoarseSpaceSolver.h"
#include "parallel/ParallelDofMapping.h"
#include "parallel/MatrixNnzStructure.h"

namespace AMDiS { namespace Parallel {

  using namespace std;

  ParallelCoarseSpaceSolver::ParallelCoarseSpaceSolver(string name)
    : ParallelSolver(name, true),
      lastMeshNnz(-1),
      alwaysCreateNnzStructure(false),
      rStartInterior(0),
      nGlobalOverallInterior(0),
      initFileStr(name)
  {
    Parameters::get("parallel->always create nnz structure", 
		    alwaysCreateNnzStructure);
    
  }


  void ParallelCoarseSpaceSolver::init(vector<const FiniteElemSpace*> &fe0,
				       vector<const FiniteElemSpace*> &fe1,
				       bool createGlobalMapping)
  {
    FUNCNAME("ParallelCoarseSpaceSolver::init()");
    
    TEST_EXIT(meshDistributor)("No mesh distributor object defined!\n");
    
    domainComm = meshDistributor->getMpiComm(meshLevel);
      
    if (meshLevel >= 1)
      coarseSpaceComm = meshDistributor->getMpiComm(meshLevel - 1);
    
    ParallelSolver::init(fe0, fe1, createGlobalMapping);
  }


  void ParallelCoarseSpaceSolver::setCoarseSpaceDofMapping(ParallelDofMapping *coarseDofs, 
							   int component)
  {
    FUNCNAME("ParallelCoarseSpaceSolver::setCoarseSpaceDofMapping()");

    TEST_EXIT_DBG(coarseDofs)("Should not happen!\n");

    if (component == -1) {
      // === Set coarse space for all components. ===

      coarseSpaceMap.clear();

      int nComponents = coarseDofs->getNumberOfComponents();
      for (int i = 0; i < nComponents; i++)
	coarseSpaceMap[i] = coarseDofs;
    } else {
      // === Set coarse space for just one component. ===

      coarseSpaceMap[component] = coarseDofs;
    }

    if (find(uniqueCoarseMap.begin(), uniqueCoarseMap.end(), coarseDofs) == 
	uniqueCoarseMap.end())
      uniqueCoarseMap.push_back(coarseDofs);
  }
  

  /// initializes vecSol, vecRhs, mat und componentIthCoarseMap
  void ParallelCoarseSpaceSolver::prepare()
  {
    FUNCNAME("ParallelCoarseSpaceSolver:prepare()");

    TEST_EXIT(uniqueCoarseMap.size() <= 2)
      ("Not yet implemented for more than two coarse spaces!\n");

    // === Create pointers to PETSc matrix and vector objects. ===

    int nCoarseMap = uniqueCoarseMap.size();
    mat.resize(nCoarseMap + 1);
    for (int i = 0; i < nCoarseMap + 1; i++)
      mat[i].resize(nCoarseMap + 1);

    vecSol.resize(nCoarseMap + 1);
    vecRhs.resize(nCoarseMap + 1);
   
    // === Create map from component number to its coarse space map. ===

    componentIthCoarseMap.resize(coarseSpaceMap.size());
    for (unsigned int i = 0; i < coarseSpaceMap.size(); i++) {
      bool found = false;
      for (int j = 0; j < nCoarseMap; j++) {
	if (coarseSpaceMap[i] == uniqueCoarseMap[j]) {
	  componentIthCoarseMap[i] = j;
	  found = true;
	  break;
	}
      }
      
      TEST_EXIT_DBG(found)("Should not happen!\n");
    }
  }


  void ParallelCoarseSpaceSolver::createMatVec(Matrix<DOFMatrix*>& seqMat)
  {
    FUNCNAME("ParallelCoarseSpaceSolver::createMatVec()");

    // === Prepare coarse space information and generate the correct number ===
    // === of empty PETSc matrix and vector objects.                        ===

    prepare();

    
    // === Update subdomain data (mostly required for multi-level methods) ===

    updateSubdomainData();


    // === If required, recompute non zero structure of the matrix. ===

    bool localMatrix = (domainComm == MPI::COMM_SELF);

    if (checkMeshChange()) {
      int nMat = uniqueCoarseMap.size() + 1; // Not multilevel should always be 1
      nnz.resize(nMat);
      for (int i = 0; i < nMat; i++) {
	nnz[i].resize(nMat);
	for (int j = 0; j < nMat; j++)
	  nnz[i][j].clear();
      }
      
      nnz[0][0].create(seqMat, *interiorMap,
		       (coarseSpaceMap.size() == 0 ? &(meshDistributor->getPeriodicMap()) : NULL),
		       meshDistributor->getElementObjectDb(),
		       localMatrix);

      for (int i = 0; i < nMat; i++) {
	for (int j = 0; j < nMat; j++) {
	  if (i == 0 && j == 0)
	    continue;
	  
	  ParallelDofMapping &rowMap = 
	    (i == 0 ? *interiorMap : *(uniqueCoarseMap[i - 1]));
	  ParallelDofMapping &colMap =
	    (j == 0 ? *interiorMap : *(uniqueCoarseMap[j - 1]));

	  nnz[i][j].create(seqMat, rowMap, colMap, NULL,
			   meshDistributor->getElementObjectDb());
	}
      }
    }


    // === Create PETSc matrices and PETSc vectors with the computed  ===
    // === nnz data structure.                                        ===
    int nRankInteriorRows = interiorMap->getRankDofs();
    int nOverallInteriorRows = interiorMap->getOverallDofs();

    if (localMatrix) {
      MatCreateSeqAIJ(domainComm, nRankInteriorRows, nRankInteriorRows,
		      0, nnz[0][0].dnnz,
		      &mat[0][0]);    
    } else {
      if (meshLevel == 0 || domainComm == MPI::COMM_SELF) {
    
	MatCreateAIJ(domainComm, nRankInteriorRows, nRankInteriorRows, 
		     nOverallInteriorRows, nOverallInteriorRows,
		     0, nnz[0][0].dnnz, 0, nnz[0][0].onnz,		   
		     &mat[0][0]);    
      } else {
	MSG("WARNING: USE MULTILEVEL, MAKE THIS GENERAL: %d of %d\n", domainComm.Get_rank(), domainComm.Get_size());
	// TODO: MAKE NNZ CALCULATION GENERAL (see also creation of coarse space
	// matrices some lines below)

	MatCreateAIJ(domainComm, nRankInteriorRows, nRankInteriorRows, 
		     nOverallInteriorRows, nOverallInteriorRows,
		     300, PETSC_NULL, 300, PETSC_NULL, &mat[0][0]);    
      }
    }
    
    MatSetOption(mat[0][0], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    if (meshLevel == 0) {
      VecCreateMPI(domainComm, nRankInteriorRows, nOverallInteriorRows, 
		   &vecSol[0]);
      VecCreateMPI(domainComm, nRankInteriorRows, nOverallInteriorRows, 
		   &vecRhs[0]);
    } else {
      MSG("WARNING: USE MULTILEVEL AND THIS IS A VERY BAD HACK 1, MAKE GENERAL!!\n");

      int nAllRows = nRankInteriorRows;
      mpi::globalAdd(coarseSpaceComm, nAllRows);
      VecCreateMPI(coarseSpaceComm, nRankInteriorRows, nAllRows, 
		   &vecSol[0]);
      VecCreateMPI(coarseSpaceComm, nRankInteriorRows, nAllRows, 
		   &vecRhs[0]);
    }

    int nCoarseMap = uniqueCoarseMap.size();
    for (int i = 0; i < nCoarseMap; i++) {
      ParallelDofMapping* cMap = uniqueCoarseMap[i];
      
      int nRankCoarseRows = cMap->getRankDofs();
      int nOverallCoarseRows = cMap->getOverallDofs();
      
      if (meshLevel == 0) {
	MatCreateAIJ(domainComm,
		     nRankCoarseRows, nRankCoarseRows,
		     nOverallCoarseRows, nOverallCoarseRows,
		     0, nnz[i + 1][i + 1].dnnz, 0, nnz[i + 1][i + 1].onnz,
		     &mat[i + 1][i + 1]);
      } else if (domainComm == MPI::COMM_SELF) {
	MatCreateAIJ(coarseSpaceComm,
		     nRankCoarseRows, nRankCoarseRows,
		     nOverallCoarseRows, nOverallCoarseRows,
		     0, nnz[i + 1][i + 1].dnnz, 0, nnz[i + 1][i + 1].onnz,
		     &mat[i + 1][i + 1]);
      } else {
	MSG("WARNING: USE MULTILEVEL AND THIS IS A VERY BAD HACK 2, MAKE GENERAL!!\n");
	
	MatCreateAIJ(coarseSpaceComm,
		     nRankCoarseRows, nRankCoarseRows,
		     nOverallCoarseRows, nOverallCoarseRows,
		     300, PETSC_NULL, 300, PETSC_NULL, 
		     &mat[i + 1][i + 1]);
      }
      MatSetOption(mat[i + 1][i + 1], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      createVec(*cMap, vecSol[i + 1]);
      createVec(*cMap, vecRhs[i + 1]);
    }

    for (int i = 0; i < nCoarseMap + 1; i++) {
      for (int j = 0; j < nCoarseMap + 1; j++) {
	if (i == j)
	  continue;

	int nRowsRankMat = (i == 0 ? nRankInteriorRows : uniqueCoarseMap[i - 1]->getRankDofs());
	int nRowsOverallMat = (i == 0 ? nOverallInteriorRows : uniqueCoarseMap[i - 1]->getOverallDofs());

	int nColsRankMat = (j == 0 ? nRankInteriorRows : uniqueCoarseMap[j - 1]->getRankDofs());
	int nColsOverallMat = (j == 0 ? nOverallInteriorRows : uniqueCoarseMap[j - 1]->getOverallDofs());

	if (meshLevel == 0) {
	  MatCreateAIJ(domainComm,
		       nRowsRankMat, nColsRankMat,
		       nRowsOverallMat, nColsOverallMat,
		       0, nnz[i][j].dnnz, 0, nnz[i][j].onnz,
		       &mat[i][j]);	  
	  MatCreateAIJ(domainComm,
		       nColsRankMat, nRowsRankMat,
		       nColsOverallMat, nRowsOverallMat,
		       0, nnz[j][i].dnnz, 0, nnz[j][i].onnz,
		       &mat[j][i]);
	} else {
	  MSG("WARNING: USE MULTILEVEL AND THIS IS A VERY BAD HACK 3, MAKE GENERAL!!\n");	  

	  if (i == 0) {
	    nRowsOverallMat = nRowsRankMat;
	    mpi::globalAdd(coarseSpaceComm, nRowsOverallMat);
	  }

	  if (j == 0) {
	    nColsOverallMat = nColsRankMat;
	    mpi::globalAdd(coarseSpaceComm, nColsOverallMat);
	  }

	  MatCreateAIJ(coarseSpaceComm,
		       nRowsRankMat, nColsRankMat,
		       nRowsOverallMat, nColsOverallMat,
		       300, PETSC_NULL, 300, PETSC_NULL,
		       &mat[i][j]);	  
	  MatCreateAIJ(coarseSpaceComm,
		       nColsRankMat, nRowsRankMat,
		       nColsOverallMat, nRowsOverallMat,
		       300, PETSC_NULL, 300, PETSC_NULL,
		       &mat[j][i]);
	}
      }
    }
  }


  void ParallelCoarseSpaceSolver::matDestroy()
  {
    FUNCNAME("ParallelCoarseSpaceSolver::matDestroy()");

    int nMatrix = mat.size();
    for (int i = 0; i < nMatrix; i++)
      for (int j = 0; j < nMatrix; j++)
	MatDestroy(&mat[i][j]);
  }


  void ParallelCoarseSpaceSolver::vecDestroy()
  {
    FUNCNAME("ParallelCoarseSpaceSolver::vecDestroy()");

    int nVec = vecSol.size();
    
    for (int i = 0; i < nVec; i++) {
      VecDestroy(&vecSol[i]);
      VecDestroy(&vecRhs[i]);
    }
  }


  void ParallelCoarseSpaceSolver::matAssembly()
  {
    FUNCNAME("ParallelCoarseSpaceSolver::matAssembly()");

    int nMatrix = mat.size();
    for (int i = 0; i < nMatrix; i++) {
      for (int j = 0; j < nMatrix; j++) {
	MatAssemblyBegin(mat[i][j], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat[i][j], MAT_FINAL_ASSEMBLY);  
      }
    }
  }


  void ParallelCoarseSpaceSolver::vecRhsAssembly()
  {
    FUNCNAME("ParallelCoarseSpaceSolver::vecRhsAssembly()");

    int nVec = vecRhs.size();
    for (int i = 0; i < nVec; i++) {
      VecAssemblyBegin(vecRhs[i]);
      VecAssemblyEnd(vecRhs[i]);
    }
  }


  void ParallelCoarseSpaceSolver::vecSolAssembly()
  {
    FUNCNAME("ParallelCoarseSpaceSolver::vecSolAssembly()");

    int nVec = vecRhs.size();
    for (int i = 0; i < nVec; i++) {
      VecAssemblyBegin(vecSol[i]);
      VecAssemblyEnd(vecSol[i]);
    }
  }


  bool ParallelCoarseSpaceSolver::checkMeshChange()
  {
    FUNCNAME("ParallelCoarseSpaceSolver::checkMeshChange()");

    int recvAllValues = 0;
    int sendValue = 
      static_cast<int>(meshDistributor->getLastMeshChangeIndex() != lastMeshNnz);
    domainComm.Allreduce(&sendValue, &recvAllValues, 1, MPI_INT, MPI_SUM);

    if (recvAllValues != 0 || alwaysCreateNnzStructure) {
      lastMeshNnz = meshDistributor->getLastMeshChangeIndex();
      return true;
    }

    return false;
  }


  void ParallelCoarseSpaceSolver::updateSubdomainData()
  {
    FUNCNAME("ParallelCoarseSpaceSolver::updateSubdomainData()");

    if (domainComm == MPI::COMM_WORLD &&
	domainComm == interiorMap->getMpiComm()) {
       rStartInterior = 0;
       nGlobalOverallInterior = interiorMap->getOverallDofs();
    } else {
      int groupRowsInterior = 0;
      if (domainComm.Get_rank() == 0)
	groupRowsInterior = interiorMap->getOverallDofs();
      
      mpi::getDofNumbering(coarseSpaceComm, groupRowsInterior,
			   rStartInterior, nGlobalOverallInterior);
      
      int tmp = 0;
      if (domainComm.Get_rank() == 0)
	tmp = rStartInterior;
      
      domainComm.Allreduce(&tmp, &rStartInterior, 1, MPI_INT, MPI_SUM);
    }
  }

} } // end namespaces
