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


#include "parallel/PetscHelper.h"
#include "parallel/PetscSolver.h"
#include "Global.h"

namespace AMDiS
{
  namespace Parallel
  {
    namespace petsc_helper
    {

      using namespace std;
      
      void getMatLocalColumn(Mat mat, PetscMatCol &matCol)
      {
	PetscInt firstIndex, lastIndex;
	MatGetOwnershipRange(mat, &firstIndex, &lastIndex);
      
	PetscInt nCols;
	const PetscInt *cols;
	const PetscScalar *values;
	
	for (int row = firstIndex; row < lastIndex; row++) {
	  MatGetRow(mat, row, &nCols, &cols, &values);
	  
	  for (int i = 0; i < nCols; i++) {
	    if (values[i] != 0.0) {
	      matCol[cols[i]].first.push_back(row - firstIndex);
	      matCol[cols[i]].second.push_back(values[i]);
	    }
	  }

	  MatRestoreRow(mat, row, &nCols, &cols, &values);
	}
      }


      void setMatLocalColumn(Mat mat, int column, Vec vec)
      {
	PetscInt firstIndex;
	MatGetOwnershipRange(mat, &firstIndex, PETSC_NULL);

	PetscInt vecSize;
	VecGetLocalSize(vec, &vecSize);

	PetscScalar *tmpValues;
	VecGetArray(vec, &tmpValues);
	for (int i  = 0; i < vecSize; i++)
	  MatSetValue(mat, 
		      firstIndex + i,
		      column,
		      tmpValues[i],
		      ADD_VALUES);
	VecRestoreArray(vec, &tmpValues);
      }


      void getColumnVec(const SparseCol &matCol, Vec vec)
      {
	VecSet(vec, 0.0);
	VecSetValues(vec, matCol.first.size(), 
		    &(matCol.first[0]), &(matCol.second[0]), INSERT_VALUES);
	VecAssemblyBegin(vec);
	VecAssemblyEnd(vec);
      }


      void blockMatMatSolve(MPI::Intracomm mpiComm, KSP ksp, Mat mat0, Mat &mat1)
      {
	FUNCNAME("blockMatMatSolve()");

	// === We have to  calculate mat1 = ksp mat0:                       ===
	// ===    - get all local column vectors from mat0                  ===
	// ===    - solve with ksp for this column vector as the rhs vector ===
	// ===    - set the result to the corresponding column vector of    ===
	// ===      matrix mat1                                             ===
	
	// Transform matrix mat0 into a sparse column format.
	PetscMatCol mat0_cols;
	getMatLocalColumn(mat0, mat0_cols);

	int nFirstCol, nLastCol;
	MatGetOwnershipRangeColumn(mat0, &nFirstCol, &nLastCol);

	int dnnz = 0;
	int onnz = 0;
	for (PetscMatCol::iterator it = mat0_cols.begin(); 
	    it != mat0_cols.end(); ++it) {
	  if (it->first >= nFirstCol && it->first < nLastCol)
	    dnnz++;
	  else
	    onnz++;
	}


	PetscInt localRows, localCols, globalRows, globalCols;
	MatGetLocalSize(mat0, &localRows, &localCols);
	MatGetSize(mat0, &globalRows, &globalCols);

	MatCreateAIJ(mpiComm,
		    localRows, localCols, globalRows, globalCols,
		    dnnz, PETSC_NULL, onnz, PETSC_NULL, &mat1);
	MatSetOption(mat1, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	
	
	Vec tmpVec;
	VecCreateSeq(PETSC_COMM_SELF, localRows, &tmpVec);

	// Solve for all column vectors of mat A_BPi and create matrix matK
	for (PetscMatCol::iterator it = mat0_cols.begin(); 
	    it != mat0_cols.end(); ++it) {
	  getColumnVec(it->second, tmpVec);
	  KSPSolve(ksp, tmpVec, tmpVec);
	  setMatLocalColumn(mat1, it->first, tmpVec);
	}
	
	VecDestroy(&tmpVec);

	MatAssemblyBegin(mat1, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat1, MAT_FINAL_ASSEMBLY);
      }


      void matNestConvert(Mat matNest, Mat &mat)  
      {
	FUNCNAME("matNestConvert()");

	PetscInt nestRows, nestCols;
	MatNestGetSize(matNest, &nestRows, &nestCols);

	TEST_EXIT(nestRows == 2 && nestCols == 2)
	  ("This funciton is only implemented for 2x2 nested matrices!\n");

	Mat mat00, mat01, mat10, mat11;
	MatNestGetSubMat(matNest, 0, 0, &mat00);
	MatNestGetSubMat(matNest, 0, 1, &mat01);
	MatNestGetSubMat(matNest, 1, 0, &mat10);
	MatNestGetSubMat(matNest, 1, 1, &mat11);

	int nRankRows0, nOverallRows0;
	MatGetLocalSize(mat00, &nRankRows0, PETSC_NULL);
	MatGetSize(mat00, &nOverallRows0, PETSC_NULL);

	int nRankRows1, nOverallRows1;
	MatGetLocalSize(mat11, &nRankRows1, PETSC_NULL);
	MatGetSize(mat11, &nOverallRows1, PETSC_NULL);

	int firstRow0;
	MatGetOwnershipRange(mat00, &firstRow0, PETSC_NULL);

	int firstRow1;
	MatGetOwnershipRange(mat11, &firstRow1, PETSC_NULL);

	int nRankRows = nRankRows0 + nRankRows1;
	int nOverallRows = nOverallRows0 + nOverallRows1;
	int firstRow = firstRow0 + firstRow1;

	int mpiSize = MPI::COMM_WORLD.Get_size();
#if (DEBUG != 0)
	int mpiRank = MPI::COMM_WORLD.Get_rank();
#endif
	vector<int> allFirstRow0(mpiSize + 1, 0);
	vector<int> allFirstRow1(mpiSize + 1, 0);
	MPI::COMM_WORLD.Allgather(&nRankRows0, 1, MPI_INT, &(allFirstRow0[1]), 1, MPI_INT);
	MPI::COMM_WORLD.Allgather(&nRankRows1, 1, MPI_INT, &(allFirstRow1[1]), 1, MPI_INT);

	for (int i = 1; i < mpiSize + 1; i++) {
	  allFirstRow0[i] += allFirstRow0[i - 1];
	  allFirstRow1[i] += allFirstRow1[i - 1];
	}

	TEST_EXIT_DBG(allFirstRow0[mpiRank] == firstRow0)("Should not happen!\n");
	TEST_EXIT_DBG(allFirstRow1[mpiRank] == firstRow1)("Should not happen!\n");

	MatCreateAIJ(PETSC_COMM_WORLD, 
		    nRankRows, nRankRows,
		    nOverallRows, nOverallRows,
		    25, PETSC_NULL, 25, PETSC_NULL, &mat);
	MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

	PetscInt nCols;
	const PetscInt *cols;
	const PetscScalar *vals;

	for (int i = 0; i < nRankRows0; i++) {
	  PetscInt newRowIndex = firstRow + i;

	  MatGetRow(mat00, firstRow0 + i, &nCols, &cols, &vals);
	  for (int j = 0; j < nCols; j++) {
	    int rank = -1;
	    for (int k = 0; k < mpiSize; k++) {
	      if (cols[j] >= allFirstRow0[k] && cols[j] < allFirstRow0[k + 1]) {
		rank = k;
		break;
	      }
	    }
	    TEST_EXIT_DBG(rank != -1)("Should not happen!\n");

	    PetscInt newColIndex = cols[j] + allFirstRow1[rank];
	    MatSetValue(mat, newRowIndex, newColIndex, vals[j], INSERT_VALUES);
	  }
	  MatRestoreRow(mat00, firstRow0 + i, &nCols, &cols, &vals);

	  MatGetRow(mat01, firstRow0 + i, &nCols, &cols, &vals);
	  for (int j = 0; j < nCols; j++) {
	    int rank = -1;
	    for (int k = 0; k < mpiSize; k++) {
	      if (cols[j] >= allFirstRow1[k] && cols[j] < allFirstRow1[k + 1]) {
		rank = k;
		break;
	      }
	    }
	    TEST_EXIT_DBG(rank != -1)("Should not happen!\n");

	    PetscInt newColIndex = cols[j] + allFirstRow0[rank + 1];
	    MatSetValue(mat, newRowIndex, newColIndex, vals[j], INSERT_VALUES);
	  }
	  MatRestoreRow(mat01, firstRow0 + i, &nCols, &cols, &vals);	
	}

	for (int i = 0; i < nRankRows1; i++) {
	  PetscInt newRowIndex = firstRow + nRankRows0 + i;

	  MatGetRow(mat10, firstRow1 + i, &nCols, &cols, &vals);
	  for (int j = 0; j < nCols; j++) {
	    int rank = -1;
	    for (int k = 0; k < mpiSize; k++) {
	      if (cols[j] >= allFirstRow0[k] && cols[j] < allFirstRow0[k + 1]) {
		rank = k;
		break;
	      }
	    }
	    TEST_EXIT_DBG(rank != -1)("Should not happen!\n");

	    PetscInt newColIndex = cols[j] + allFirstRow1[rank];
	    MatSetValue(mat, newRowIndex, newColIndex, vals[j], INSERT_VALUES);
	  }
	  MatRestoreRow(mat10, firstRow1 + i, &nCols, &cols, &vals);

	  MatGetRow(mat11, firstRow1 + i, &nCols, &cols, &vals);
	  for (int j = 0; j < nCols; j++) {
	    int rank = -1;
	    for (int k = 0; k < mpiSize; k++) {
	      if (cols[j] >= allFirstRow1[k] && cols[j] < allFirstRow1[k + 1]) {
		rank = k;
		break;
	      }
	    }
	    TEST_EXIT_DBG(rank != -1)("Should not happen!\n");

	    PetscInt newColIndex = cols[j] + allFirstRow0[rank + 1];
	    MatSetValue(mat, newRowIndex, newColIndex, vals[j], INSERT_VALUES);
	  }
	  MatRestoreRow(mat11, firstRow1 + i, &nCols, &cols, &vals);
	}

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
      }


      void setSolverWithLu(KSP ksp, 
			  const char* kspPrefix,
			  KSPType kspType, 
			  PCType pcType, 
			  const MatSolverPackage matSolverPackage,
			  PetscReal rtol,
			  PetscReal atol,
			  PetscInt maxIt)
      {
	KSPSetType(ksp, kspType);
	KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT, maxIt);
	//if (kspPrefix != "")
	if (strcmp(kspPrefix, "") != 0)
	  KSPSetOptionsPrefix(ksp, kspPrefix);
	KSPSetFromOptions(ksp);

	PC pc;
	KSPGetPC(ksp, &pc);
	PCSetType(pc, pcType);
	if (matSolverPackage != PETSC_NULL)
	  PCFactorSetMatSolverPackage(pc, matSolverPackage);
	PCSetFromOptions(pc);
	
#if DEBUG != 0
    MSG("PetscOptionsView:\n");
    PetscViewer viewer;
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERASCII);
    PetscOptionsView(viewer);
    PetscViewerDestroy(&viewer);
#endif
      }


      void setSolver(KSP ksp, 
		    const char* kspPrefix,
		    KSPType kspType, 
		    PCType pcType, 
		    PetscReal rtol,
		    PetscReal atol,
		    PetscInt maxIt)
      {

	setSolverWithLu(ksp, kspPrefix, kspType, pcType, PETSC_NULL, 
			rtol, atol, maxIt);
      }
      
      
      void createSolver(MPI::Intracomm comm, KSP &ksp, Mat m, std::string kspPrefix, int info)
      {
	KSPCreate(comm, &ksp);
	#if (PETSC_VERSION_MINOR >= 5)
	  KSPSetOperators(ksp, m, m);
	#else
	  KSPSetOperators(ksp, m, m, SAME_NONZERO_PATTERN);
	#endif
	KSPSetTolerances(ksp, 0.0, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetType(ksp, KSPBCGS);
	KSPSetOptionsPrefix(ksp, kspPrefix.c_str());
	
	if (info >= 10)
	  KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL);
	else if (info >= 20)
	  KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, PETSC_NULL, PETSC_NULL);
      }
      
    } // end namespace petsc_helper
  } // end namespace Parallel
} // end namespace AMDiS
