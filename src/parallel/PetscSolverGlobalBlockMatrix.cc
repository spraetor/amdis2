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


#include "parallel/PetscSolverGlobalBlockMatrix.h"
#include "parallel/StdMpi.h"
#include "parallel/MpiHelper.h"

using namespace std;

namespace AMDiS { namespace Parallel {

  void PetscSolverGlobalBlockMatrix::fillPetscMatrix(Matrix<DOFMatrix*> *seqMat)
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::fillPetscMatrix()");

    TEST_EXIT_DBG(meshDistributor)("No mesh distributor object defined!\n");
    TEST_EXIT_DBG(interiorMap)("No parallel mapping object defined!\n");
    TEST_EXIT_DBG(seqMat)("No DOF matrix defined!\n");

    double wtime = MPI::Wtime();

    prepare();

    const FiniteElemSpace *feSpace = componentSpaces[0];
    nComponents = seqMat->getNumRows();
    int nRankRows = (*interiorMap)[feSpace].nRankDofs;
    int nOverallRows = (*interiorMap)[feSpace].nOverallDofs;

#if (DEBUG != 0)
    MSG("Fill petsc matrix 1 needed %.5f seconds\n", MPI::Wtime() - wtime);
#endif

    if (nBlocks == -1) {
      nBlocks = nComponents;
      for (int i = 0; i < nBlocks; i++)
	componentInBlock[i] = i;
    }

    vector<int> compNthInBlock(nComponents, 0);
    vector<int> blockSize(nBlocks, 0);

    for (int i = 0; i < nComponents; i++) {
      compNthInBlock[i] = blockSize[componentInBlock[i]];
      blockSize[componentInBlock[i]]++;
    }

    nestMat.resize(nBlocks * nBlocks);

    // === Transfer values from DOF matrices to the PETSc matrix. === 

    for (int i = 0; i < nBlocks; i++)
      for (int j = 0; j < nBlocks; j++)
	MatCreateAIJ(domainComm,
		     nRankRows * blockSize[i], nRankRows * blockSize[j],
		     nOverallRows * blockSize[i], nOverallRows * blockSize[j],
		     300 * blockSize[i], PETSC_NULL, 
		     300 * blockSize[j], PETSC_NULL,
		     &(nestMat[i * nBlocks + j]));
			
    for (int i = 0; i < nComponents; i++)
      for (int j = 0; j < nComponents; j++)
	if ((*seqMat)[i][j]) {
	  int idx = componentInBlock[i] * nBlocks + componentInBlock[j];
	  setDofMatrix(nestMat[idx], (*seqMat)[i][j], 
		       compNthInBlock[i], compNthInBlock[j]);
	}

    for (int i = 0; i < nBlocks; i++) {
      for (int j = 0; j < nBlocks; j++) {
	int idx = i * nBlocks + j;
	if (nestMat[idx]) {
	  MatAssemblyBegin(nestMat[idx], MAT_FINAL_ASSEMBLY);
	  MatAssemblyEnd(nestMat[idx], MAT_FINAL_ASSEMBLY);
	}
      }
    }	  
	

    MatCreateNest(domainComm, nBlocks, PETSC_NULL, nBlocks, PETSC_NULL,
		  &(nestMat[0]), &getMatInterior());

#if (DEBUG != 0)
    MSG("Fill petsc matrix 2 needed %.5f seconds\n", MPI::Wtime() - wtime);
#endif

    matAssembly();

    
    /// initPreconditioner(...)

    // === Init PETSc solver and preconditioner objects. ===

    initSolver(kspInterior);
    KSPGetPC(kspInterior, &pcInterior);
    initPreconditioner(pcInterior);

    MSG("Fill petsc matrix needed %.5f seconds\n", MPI::Wtime() - wtime);
  }


  void PetscSolverGlobalBlockMatrix::fillPetscRhs(SystemVector *vec)
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::fillPetscRhs()");

    TEST_EXIT_DBG(vec)("NO DOF vector defined!\n");

    nComponents = vec->getSize();
    const FiniteElemSpace *feSpace = componentSpaces[0];
    int nRankRows = (*interiorMap)[feSpace].nRankDofs;
    int nOverallRows = (*interiorMap)[feSpace].nOverallDofs;

    nestVec.resize(nComponents);

    for (int i = 0; i < nComponents; i++) {
      VecCreateMPI(domainComm, nRankRows, nOverallRows, &(nestVec[i]));

      setDofVector(nestVec[i], vec->getDOFVector(i));
      
      VecAssemblyBegin(nestVec[i]);
      VecAssemblyEnd(nestVec[i]);
    }

    VecCreateNest(domainComm, nComponents, PETSC_NULL, 
		  &(nestVec[0]), &(getVecRhsInterior()));

    vecRhsAssembly();
  }


  void PetscSolverGlobalBlockMatrix::initSolver(KSP &ksp)
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::initSolver()");

    KSPCreate(domainComm, &ksp);
#if (PETSC_VERSION_MINOR >= 5)
    KSPSetOperators(ksp, getMatInterior(), getMatInterior());
#else
    KSPSetOperators(ksp, getMatInterior(), getMatInterior(), SAME_NONZERO_PATTERN);
#endif
    KSPSetOptionsPrefix(ksp, kspPrefix.c_str());
    KSPSetFromOptions(ksp);
  }


  void PetscSolverGlobalBlockMatrix::exitSolver(KSP ksp)
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::exitSolver()");

    KSPDestroy(&ksp);
  }


  void PetscSolverGlobalBlockMatrix::initPreconditioner(PC pc)
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::initPreconditioner()");

    PCSetFromOptions(pc);
  }


  void PetscSolverGlobalBlockMatrix::exitPreconditioner(PC pc)
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::exitPreconditioner()");
  }


  void PetscSolverGlobalBlockMatrix::solvePetscMatrix(SystemVector &vec, 
						      AdaptInfo *adaptInfo)
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::solvePetscMatrix()");

    const FiniteElemSpace *feSpace = componentSpaces[0];
    VecDuplicate(getVecRhsInterior(), &petscSolVec);
    
    for (int i = 0; i < vec.getSize(); i++)
    {
      Vec tmp;
      VecNestGetSubVec(petscSolVec, i, &tmp);
      setDofVector(tmp, vec.getDOFVector(i));
      VecAssemblyBegin(tmp); 
      VecAssemblyEnd(tmp);
    }

    // PETSc.
    solve(getVecRhsInterior(), petscSolVec);

    // === Transfere values from PETSc's solution vectors to the DOF vectors. ===
    for (int i = 0; i < nComponents; i++) {
      DOFVector<double> &dofvec = *(vec.getDOFVector(i));

      Vec tmp;
      VecNestGetSubVec(petscSolVec, i, &tmp);

//       int nRankDofs = (*interiorMap)[feSpace].nRankDofs;
      PetscScalar *vecPointer;
      VecGetArray(tmp, &vecPointer);

      DofMap& d = (*interiorMap)[feSpace].getMap();
      for (DofMap::iterator it = d.begin(); it != d.end(); ++it)
	if (it->second.local != -1)
	  dofvec[it->first] = vecPointer[it->second.local];

      VecRestoreArray(tmp, &vecPointer);
    }

    VecDestroy(&petscSolVec);

    // === Synchronize DOFs at common DOFs, i.e., DOFs that correspond to ===
    // === more than one partition.                                       ===
    meshDistributor->synchVector(vec);
  }


  void PetscSolverGlobalBlockMatrix::destroyMatrixData()
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::destroyMatrixData()");

    for (unsigned int i = 0; i < nestMat.size(); i++)
      if (nestMat[i] != PETSC_NULL)
	MatDestroy(&(nestMat[i]));

    matDestroy();

    exitPreconditioner(pcInterior);

    exitSolver(kspInterior);
  }


  void PetscSolverGlobalBlockMatrix::destroyVectorData()
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::destroyVectorData()");

    for (unsigned int i = 0; i < nestVec.size(); i++)
      VecDestroy(&(nestVec[i]));
      
    vecDestroy();
  }


  void PetscSolverGlobalBlockMatrix::setDofMatrix(Mat& petscMat, 
						  DOFMatrix* seqMat,
						  int dispRowBlock, 
						  int dispColBlock)
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::setDofMatrix()");

    TEST_EXIT(petscMat)("No PETSc matrix!\n");
    TEST_EXIT(seqMat)("No DOFMatrix!\n");

    const FiniteElemSpace *feSpace = componentSpaces[0];

    using mtl::tag::row; using mtl::tag::nz; using mtl::begin; using mtl::end;
    namespace traits = mtl::traits;
    typedef DOFMatrix::base_matrix_type Matrix;

    traits::col<Matrix>::type col(seqMat->getBaseMatrix());
    traits::const_value<Matrix>::type value(seqMat->getBaseMatrix());

    typedef traits::range_generator<row, Matrix>::type cursor_type;
    typedef traits::range_generator<nz, cursor_type>::type icursor_type;

    int dispRowIndex = (*interiorMap)[feSpace].nRankDofs * dispRowBlock;
    int dispColIndex = (*interiorMap)[feSpace].nRankDofs * dispColBlock;

    vector<int> cols;
    vector<double> values;
    cols.reserve(3000);
    values.reserve(3000);
    
    // === Traverse all rows of the dof matrix and insert row wise the values ===
    // === to the PETSc matrix.                                               ===

    for (cursor_type cursor = begin<row>(seqMat->getBaseMatrix()), 
	   cend = end<row>(seqMat->getBaseMatrix()); cursor != cend; ++cursor) {

      // Global index of the current row DOF.
      int rowIndex = (*interiorMap)[feSpace][cursor.value()].global + dispRowIndex;

      cols.clear();
      values.clear();
      
      for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor); 
	   icursor != icend; ++icursor) {	
	// Global index of the current column index.
	int colIndex = (*interiorMap)[feSpace][col(*icursor)].global + dispColIndex;
	
	// Ignore all zero entries, expect it is a diagonal entry.
	if (value(*icursor) == 0.0 && rowIndex != colIndex)
	  continue;
	
	// Calculate the exact position of the column index in the PETSc matrix.
	cols.push_back(colIndex);
	values.push_back(value(*icursor));
      }

      MatSetValues(petscMat, 1, &rowIndex, cols.size(), 
 		   &(cols[0]), &(values[0]), ADD_VALUES);	
    }
  }
  
  
  void PetscSolverGlobalBlockMatrix::setDofVector(Vec& petscVec, 
						  DOFVector<double>* vec)
  {
    FUNCNAME("PetscSolverGlobalBlockMatrix::setDofVector()");

    const FiniteElemSpace *feSpace = componentSpaces[0];

    // Traverse all used DOFs in the dof vector.
    DOFVector<double>::Iterator dofIt(vec, USED_DOFS);
    for (dofIt.reset(); !dofIt.end(); ++dofIt) {
      int index = (*interiorMap)[feSpace][dofIt.getDOFIndex()].global;
      double value = *dofIt;

      VecSetValues(petscVec, 1, &index, &value, ADD_VALUES);
    }
  }

} } // end namespace Parallel, AMDiS
