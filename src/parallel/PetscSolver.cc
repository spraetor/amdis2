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


#include "AMDiS.h"
#include "DOFMatrix.h"
#include "Global.h"
#include "Initfile.h"
#include "parallel/MeshDistributor.h"
#include "parallel/MpiHelper.h"
#include "parallel/ParallelDofMapping.h"
#include "parallel/PetscSolver.h"
#include "parallel/StdMpi.h"

namespace AMDiS { namespace Parallel {

  PetscSolver::PetscSolver(std::string name)
    : super(name),
      kspPrefix(""),
      removeRhsNullspace(false),    
      hasConstantNullspace(false),
      isSymmetric(false),
      handleDirichletRows(true)
  {
    std::string tmp("");
    Parameters::get(name + "->petsc prefix", tmp);
    setKspPrefix(tmp); 
    
    std::string kspStr = "";
    Parameters::get("parallel->solver->petsc->ksp", kspStr);
    if (!kspStr.size())
      Parameters::get(name + "->ksp", kspStr);
      
    if (kspStr != "")
      PetscOptionsInsertString(kspStr.c_str());

    Parameters::get(name + "->remove rhs null space", removeRhsNullspace);
    Parameters::get(name + "->has constant null space", hasConstantNullspace);
    Parameters::get(name + "->nullspace->const in comp", 
		    constNullspaceComponent);
  }


  PetscSolver::~PetscSolver()
  {}


  void PetscSolver::init(std::vector<const FiniteElemSpace*> &fe0,
			 std::vector<const FiniteElemSpace*> &fe1,
			 bool createGlobalMapping)
  {
    FUNCNAME("PetscSolver::init()");

    TEST_EXIT(meshDistributor)("No mesh distributor object defined!\n");
    meshDistributor->setBoundaryDofRequirement(getBoundaryDofRequirement());

    super::init(fe0, fe1, createGlobalMapping);
  }

  /// see seq::LinearSolver::solveLinearSystem(...)
  int PetscSolver::solveLinearSystem(const SolverMatrix<Matrix<DOFMatrix*> >& A,
				    SystemVector& x,
				    SystemVector& b,
				    bool createMatrixData,
				    bool storeMatrixData)
  {
    FUNCNAME("PetscSolver::solveLinearSystem()");
    
    TEST_EXIT(meshDistributor)("No meshDistributor provided. Should not happen!\n");

    MPI::COMM_WORLD.Barrier();
    Timer t;
    
    if (createMatrixData)
      fillPetscMatrix(const_cast< Matrix< DOFMatrix* >* >(A.getOriginalMat()));    
    
    fillPetscRhs(&b);

    INFO(info, 8)("creation of parallel data structures needed %.5f seconds\n", 
		  t.elapsed());

    solvePetscMatrix(x, NULL);

    if (!storeMatrixData) {
      destroyVectorData();
      destroyMatrixData();
    }
   
    return getErrorCode();
  }

  void PetscSolver::fillPetscMatrix(DOFMatrix* mat)
  {
    Matrix<DOFMatrix*> sysMat(1, 1);
    sysMat[0][0] = mat;
    fillPetscMatrix(&sysMat);
  }


  /// see seq::PetscSolver::solve(const MatrixType& A, VectorType& x, const VectorType& b) 
  void PetscSolver::solve(Vec &rhs, Vec &sol)
  {
    PetscErrorCode solverError = KSPSolve(kspInterior, rhs, sol);
    
    PetscInt nIter = 0;
    KSPGetIterationNumber(kspInterior, &nIter);
    PetscReal residual_norm = -1.0;
    KSPGetResidualNorm(kspInterior, &residual_norm);
    
    if (residual_norm <= 0.0) {
      Vec r;
      createVec(*interiorMap, r);
      KSPBuildResidual(kspInterior, PETSC_NULL, PETSC_NULL, &r);
      VecNorm(r, NORM_2, &residual_norm);
    }
    
    setErrorCode(solverError);
    setIterations(nIter);
    setResidual(residual_norm);
  }


  /// was soll das denn genau machen?
  void PetscSolver::solveGlobal(Vec &rhs, Vec &sol)
  {
    FUNCNAME("PetscSolver::solveGlobal()");

    int s, ls;
    VecGetSize(rhs, &s);
    VecGetLocalSize(rhs, &ls);

    MSG("Solve global %d %d\n", ls, s);

    ERROR_EXIT("Not implemented!\n");
  }


  /// -> Helper-Funktion, benötigt nur domainComm
  void PetscSolver::copyVec(Vec& originVec, Vec& destVec, 
			    std::vector<int>& originIndex, std::vector<int>& destIndex)
  {
    IS originIs, destIs;
    ISCreateGeneral(domainComm, 
		    originIndex.size(), 
		    &(originIndex[0]),
		    PETSC_USE_POINTER,
		    &originIs);

    ISCreateGeneral(domainComm, 
		    destIndex.size(), 
		    &(destIndex[0]),
		    PETSC_USE_POINTER,
		    &destIs);

    VecScatter scatter;
    VecScatterCreate(originVec, originIs, destVec, destIs, &scatter);
    VecScatterBegin(scatter, originVec, destVec,
		    INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter, originVec, destVec,
		  INSERT_VALUES, SCATTER_FORWARD);

    ISDestroy(&originIs);
    ISDestroy(&destIs);    
    VecScatterDestroy(&scatter);
  }


  /// -> Helper-Funktion
  bool PetscSolver::testMatrixSymmetric(Mat mat, bool advancedTest)
  {
    FUNCNAME("PetscSolver::testMatrixSymmetric()");

    Mat matTrans;
    MatTranspose(mat, MAT_INITIAL_MATRIX, &matTrans);

    bool isSym = true;

    if (advancedTest) {
      int rowStart, rowEnd;
      MatGetOwnershipRange(mat, &rowStart, &rowEnd);

      PetscInt nCols, nColsTrans;
      const PetscInt *cols, *colsTrans;
      const PetscScalar *vals, *valsTrans;

      for (int i = rowStart; i < rowEnd; i++) {
	MatGetRow(mat, i, &nCols, &cols, &vals);
	MatGetRow(matTrans, i, &nColsTrans, &colsTrans, &valsTrans);

	if (nCols != nColsTrans) {
	  MSG("Symmetric test fails: mat row %d: %d nnz   mat col %d: %d nnz\n",
	      i, nCols, i, nColsTrans);
	  isSym = false;
	} else {
	  for (int j = 0; j < nCols; j++) {
	    if (cols[j] != colsTrans[j]) {
	      if (cols[j] > colsTrans[j]) {
		MSG("mat[%d][%d] does not exists: mat[%d][%d] = %e\n", 
		    i, colsTrans[j],
		    colsTrans[j], i, valsTrans[j]);
	      } else {
		MSG("mat[%d][%d] does not exists: mat[%d][%d] = %e\n", 
		    cols[j], i,
		    i, cols[j], vals[j]);
	      }
	      isSym = false;
	    } else {
	      double n = fabs(vals[j] - valsTrans[j]);
	      if (n > 1e-10) {
		MSG("value diff:  mat[%d][%d] = %e   mat[%d][%d] = %e\n",
		    i, cols[j], vals[j], colsTrans[j], i, valsTrans[j]);

		isSym = false;
	      }
	    }
	  }
	}

	MatRestoreRow(mat, i, &nCols, &cols, &vals);
	MatRestoreRow(matTrans, i, &nColsTrans, &colsTrans, &valsTrans);
      }
    } 
      
    MatAXPY(matTrans, -1.0, mat, DIFFERENT_NONZERO_PATTERN);
    double norm1, norm2;
    MatNorm(matTrans, NORM_FROBENIUS, &norm1);
    MatNorm(matTrans, NORM_INFINITY, &norm2);
    MatDestroy(&matTrans);
    
    MSG("Matrix norm test: %e %e\n", norm1, norm2);
    
    return (isSym && norm1 < 1e-10 && norm2 < 1e-10);
  }
} } // end namespaces
