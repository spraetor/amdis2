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
// #include "DirichletBC.h"
#include "DOFVector.h"
#include "parallel/PetscSolverGlobalMatrix.hpp"
#include "parallel/StdMpi.hpp"
#include "parallel/MpiHelper.hpp"
#include "solver/PetscTypes.hpp"

using namespace std;

namespace AMDiS
{
  namespace Parallel
  {

    PetscSolverGlobalMatrix::PetscSolverGlobalMatrix(string name, bool setOptions)
      : PetscSolver(name),
        zeroStartVector(false),
        printMatInfo(false)
    {
      FUNCNAME_DBG("PetscSolverGlobalMatrix()");

      bool matSolverPackage = false;
      if (setOptions)
      {
        PetscParameters params;

        // set the solver
        std::string solverName = "petsc";
        Parameters::get(name, solverName);
        if (solverName == "petsc")
          Parameters::get(name + "->ksp_type", solverName);

        std::string kspSolver = params.solverMap[solverName];

        if (params.matSolverPackage.find(kspSolver) != params.matSolverPackage.end())
        {
          // direct solvers
          PetscOptionsInsertString(("-" + kspPrefix + "ksp_type preonly").c_str());
          PetscOptionsInsertString(("-" + kspPrefix + "pc_type lu").c_str());
          PetscOptionsInsertString(("-" + kspPrefix + "pc_factor_mat_solver_package " + kspSolver).c_str());
          setMaxIterations(1);
          zeroStartVector = true;
          matSolverPackage = true;
        }
        else if (params.emptyParam.find(kspSolver) == params.emptyParam.end() && solverName != "petsc")
        {
          // other solvers
          PetscOptionsInsertString(("-" + kspPrefix + "ksp_type " + kspSolver).c_str());
        }

        // set the preconditioner
        string precon = "";
        Parameters::get(name + "->pc_type", precon);
        if (!precon.size())
          Parameters::get(name + "->left precon", precon);
        if (!precon.size())
          Parameters::get(name + "->right precon", precon);
        if (!matSolverPackage && params.emptyParam.find(precon) == params.emptyParam.end())
        {
          precon = (params.preconMap.find(precon) != params.preconMap.end() ? params.preconMap[precon] : precon);
          PetscOptionsInsertString(("-" + kspPrefix + "pc_type " + precon).c_str());
        }

        PetscOptionsInsertString(("-" + kspPrefix + "ksp_max_it " + std::to_string(getMaxIterations())).c_str());
        PetscOptionsInsertString(("-" + kspPrefix + "ksp_rtol " + std::to_string(getRelative())).c_str());
        PetscOptionsInsertString(("-" + kspPrefix + "ksp_atol " + std::to_string(getTolerance())).c_str());

        if (getInfo() >= 20)
          PetscOptionsInsertString(("-" + kspPrefix + "ksp_monitor_true_residual").c_str());
        else if (getInfo() >= 10)
          PetscOptionsInsertString(("-" + kspPrefix + "ksp_monitor").c_str());
      }
      if (!matSolverPackage)
      {
        Parameters::get(name + "->use zero start vector", zeroStartVector);
      }
      Parameters::get("parallel->print matrix info", printMatInfo);

#if DEBUG != 0
      bool printOptionsInfo = false;
      Parameters::get("parallel->debug->print options info", printOptionsInfo);
      if (printOptionsInfo)
      {
        MSG("PetscOptionsView:\n");
        PetscViewer viewer;
        PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
        PetscViewerSetType(viewer, PETSCVIEWERASCII);
        PetscOptionsView(viewer);
        PetscViewerDestroy(&viewer);
      }
#endif
    }


    void PetscSolverGlobalMatrix::fillPetscMatrix(Matrix<DOFMatrix*>* seqMat)
    {
      FUNCNAME("PetscSolverGlobalMatrix::fillPetscMatrix()");

      TEST_EXIT_DBG(meshDistributor)("No mesh distributor object defined!\n");
      TEST_EXIT_DBG(interiorMap)("No parallel mapping object defined!\n");
      TEST_EXIT_DBG(seqMat)("No DOF matrix defined!\n");

#if (DEBUG != 0)
      Timer t;
#endif

      createMatVec(*seqMat);

      if (coarseSpaceMap.size())
      {
        fillPetscMatrixWithCoarseSpace(seqMat);
        return;
      }

      // === Create PETSc vector (solution and a temporary vector). ===

#if (DEBUG != 0)
      MSG("Fill petsc matrix 1 needed %.5f seconds\n", t.elapsed());
      t.reset();
#endif

      // === Transfer values from DOF matrices to the PETSc matrix. ===

      int nComponents = seqMat->getNumRows();
      for (int i = 0; i < nComponents; i++)
        for (int j = 0; j < nComponents; j++)
          if ((*seqMat)[i][j])
            setDofMatrix((*seqMat)[i][j], i, j);

#if (DEBUG != 0)
      MSG("Fill petsc matrix 2 needed %.5f seconds\n", t.elapsed());
      t.reset();
#endif

      matAssembly();

      if (printMatInfo)
      {
        MatInfo matInfo;
        MatGetInfo(getMatInterior(), MAT_GLOBAL_SUM, &matInfo);
        MSG("Matrix info:\n");
        MSG("  memory usage: %e MB\n", matInfo.memory / (1024.0 * 1024.0));
        MSG("  mallocs: %d\n", static_cast<int>(matInfo.mallocs));
        MSG("  nz allocated: %d\n", static_cast<int>(matInfo.nz_allocated));
        MSG("  nz used: %d\n", static_cast<int>(matInfo.nz_used));
        MSG("  nz unneeded: %d\n", static_cast<int>(matInfo.nz_unneeded));
      }


      // === Init PETSc solver and preconditioner objects. ===

      initSolver(kspInterior);
      KSPGetPC(kspInterior, &pcInterior);
      initPreconditioner(*seqMat, mat[0][0]);


#if (DEBUG != 0)
      MSG("Fill petsc matrix 3 needed %.5f seconds\n", t.elapsed());
#endif


      // === For debugging allow to write the matrix to a file. ===

      bool dbgWriteMatrix = false;
      Parameters::get("parallel->debug->write matrix", dbgWriteMatrix);
      if (dbgWriteMatrix)
      {
        PetscViewer matView;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mpi.mat",
                              FILE_MODE_WRITE, &matView);
        MatView(getMatInterior(), matView);
        PetscViewerDestroy(&matView);
      }
    }


    void PetscSolverGlobalMatrix::fillPetscMatrixWithCoarseSpace(Matrix<DOFMatrix*>* seqMat)
    {
      FUNCNAME_DBG("PetscSolverGlobalMatrix::fillPetscMatrixWithCoarseSpace()");

      TEST_EXIT_DBG(interiorMap)("No interiorMap! Should not happen!\n");
      TEST_EXIT_DBG(coarseSpaceMap.size() == (size_t)seqMat->getSize())
      ("Wrong sizes %d %d\n", coarseSpaceMap.size(), seqMat->getSize());

      // === Prepare traverse of sequentially created matrices. ===

      using mtl::tag::row;
      using mtl::tag::nz;
      using mtl::begin;
      using mtl::end;
      namespace traits = mtl::traits;
      typedef DOFMatrix::base_matrix_type Matrix;

      typedef traits::range_generator<row, Matrix>::type cursor_type;
      typedef traits::range_generator<nz, cursor_type>::type icursor_type;

      vector<int> cols, colsOther;
      vector<double> values, valuesOther;
      cols.reserve(300);
      colsOther.reserve(300);
      values.reserve(300);
      valuesOther.reserve(300);

      bool localMatrix =
        (meshDistributor->getMeshLevelData().getMpiComm(meshLevel) == MPI::COMM_SELF);

      // === Traverse all sequentially created matrices and add the values to ===
      // === the global PETSc matrices.                                       ===

      int nComponents = seqMat->getSize();
      for (int rowComponent = 0; rowComponent < nComponents; rowComponent++)
      {
        for (int colComponent = 0; colComponent < nComponents; colComponent++)
        {
          DOFMatrix* dofMat = (*seqMat)[rowComponent][colComponent];

          if (!dofMat)
            continue;

          ParallelDofMapping* rowCoarseSpace = coarseSpaceMap[rowComponent];
          ParallelDofMapping* colCoarseSpace = coarseSpaceMap[colComponent];

          std::set<DegreeOfFreedom>& dirichletRows = dofMat->getDirichletRows();

          traits::col<Matrix>::type col(dofMat->getBaseMatrix());
          traits::const_value<Matrix>::type value(dofMat->getBaseMatrix());

          // Traverse all rows.
          for (cursor_type cursor = begin<row>(dofMat->getBaseMatrix()),
               cend = end<row>(dofMat->getBaseMatrix()); cursor != cend; ++cursor)
          {

            bool isRowCoarse = isCoarseSpace(rowComponent, cursor.value());

            // For the case, this is a dirichlet row we have to check whether the
            // rank is also owner of this row DOF.
            if (dirichletRows.count(cursor.value()))
            {
              if ((!isRowCoarse && !(*interiorMap)[rowComponent].isRankDof(cursor.value())) ||
                  (isRowCoarse && !(*rowCoarseSpace)[rowComponent].isRankDof(cursor.value())))
                continue;
            }

            cols.clear();
            colsOther.clear();
            values.clear();
            valuesOther.clear();

            // Traverse all columns.
            for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor);
                 icursor != icend; ++icursor)
            {

              bool isColCoarse = isCoarseSpace(colComponent, col(*icursor));

              if (isColCoarse == false)
                if ((*interiorMap)[colComponent].isSet(col(*icursor)) == false)
                  continue;

              if (isColCoarse == isRowCoarse)
              {
                cols.push_back(col(*icursor));
                values.push_back(value(*icursor));
              }
              else
              {
                colsOther.push_back(col(*icursor));
                valuesOther.push_back(value(*icursor));
              }
            }  // for each nnz in row


            // === Set matrix values. ===

            if (isRowCoarse)
            {
              for (unsigned int i = 0; i < cols.size(); i++)
                cols[i] = colCoarseSpace->getMatIndex(colComponent, cols[i]);

              int rowIndex = rowCoarseSpace->getMatIndex(rowComponent, cursor.value());
              MatSetValues(getMatCoarseByComponent(rowComponent, colComponent),
                           1, &rowIndex, cols.size(),
                           &(cols[0]), &(values[0]), ADD_VALUES);

              if (colsOther.size())
              {
                for (unsigned int i = 0; i < colsOther.size(); i++)
                  colsOther[i] =
                    interiorMap->getMatIndex(colComponent, colsOther[i]) +
                    rStartInterior;

                MatSetValues(getMatCoarseInteriorByComponent(rowComponent),
                             1, &rowIndex, colsOther.size(),
                             &(colsOther[0]), &(valuesOther[0]), ADD_VALUES);
              }
            }
            else
            {
              if ((*interiorMap)[rowComponent].isSet(cursor.value()) == false)
                continue;

              int localRowIndex =
                (localMatrix ?
                 interiorMap->getLocalMatIndex(rowComponent, cursor.value()) :
                 interiorMap->getMatIndex(rowComponent, cursor.value()));

              for (unsigned int i = 0; i < cols.size(); i++)
              {
                if (localMatrix)
                  cols[i] = interiorMap->getLocalMatIndex(colComponent, cols[i]);
                else
                  cols[i] = interiorMap->getMatIndex(colComponent, cols[i]);
              }

              MatSetValues(getMatInterior(), 1, &localRowIndex, cols.size(),
                           &(cols[0]), &(values[0]), ADD_VALUES);

              if (colsOther.size())
              {
                int globalRowIndex =
                  interiorMap->getMatIndex(rowComponent, cursor.value()) + rStartInterior;

                for (unsigned int i = 0; i < colsOther.size(); i++)
                  colsOther[i] =
                    colCoarseSpace->getMatIndex(colComponent, colsOther[i]);

                MatSetValues(getMatInteriorCoarseByComponent(colComponent),
                             1, &globalRowIndex, colsOther.size(),
                             &(colsOther[0]), &(valuesOther[0]), ADD_VALUES);
              }
            }
          }
        }
      }

      matAssembly();

      // === Create solver for the non primal (thus local) variables. ===

      KSPCreate(domainComm, &kspInterior);
#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(kspInterior, getMatInterior(), getMatInterior());
#else
      KSPSetOperators(kspInterior, getMatInterior(), getMatInterior(), SAME_NONZERO_PATTERN);
#endif
      KSPSetOptionsPrefix(kspInterior, "interior_");
      KSPSetType(kspInterior, KSPPREONLY);
      KSPGetPC(kspInterior, &pcInterior);
      if (isSymmetric)
      {
        PCSetType(pcInterior, PCCHOLESKY);
        PCFactorSetMatSolverPackage(pcInterior, MATSOLVERMUMPS);
      }
      else
      {
        PCSetType(pcInterior, PCLU);
        if (localMatrix)
          PCFactorSetMatSolverPackage(pcInterior, MATSOLVERUMFPACK);
        else
          PCFactorSetMatSolverPackage(pcInterior, MATSOLVERMUMPS);
      }
      KSPSetFromOptions(kspInterior);
    }


    void PetscSolverGlobalMatrix::fillPetscRhs(SystemVector* vec)
    {
      FUNCNAME_DBG("PetscSolverGlobalMatrix::fillPetscRhs()");

      TEST_EXIT_DBG(vec)("No DOF vector defined!\n");
      TEST_EXIT_DBG(interiorMap)("No parallel DOF map defined!\n");

      // === Transfer values from DOF vector to the PETSc vector. ===
      if (coarseSpaceMap.size())
      {
        for (int i = 0; i < vec->getSize(); i++)
          setDofVector(getVecRhsInterior(),
                       getVecRhsCoarseByComponent(i), vec->getDOFVector(i), i);
      }
      else
      {
        for (int i = 0; i < vec->getSize(); i++)
          setDofVector(getVecRhsInterior(), vec->getDOFVector(i), i);
      }

      vecRhsAssembly();

      // === For debugging allow to write the rhs vector to a file. ===

      bool dbgWriteRhs = false;
      Parameters::get("parallel->debug->write rhs", dbgWriteRhs);
      if (dbgWriteRhs)
      {
        PetscViewer vecView;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mpi.vec",
                              FILE_MODE_WRITE, &vecView);
        VecView(getVecRhsInterior(), vecView);
        PetscViewerDestroy(&vecView);
      }
    }


    /// 1.) set startsolution
    /// 2.) create null-space
    /// 3.) solve Ax=b
    /// 4.) destroy null-space
    /// 5.) transfer solution back to DOFVector
    void PetscSolverGlobalMatrix::solvePetscMatrix(SystemVector& vec,
        AdaptInfo& adaptInfo)
    {
      FUNCNAME("PetscSolverGlobalMatrix::solvePetscMatrix()");

      int nComponents = vec.getSize();

      // === Set old solution to be initiual guess for PETSc solver. ===
      if (!zeroStartVector)
      {
        TEST_EXIT(coarseSpaceMap.size() == 0)("Not yet supported!\n");

        VecSet(getVecSolInterior(), 0.0);

        for (int i = 0; i < nComponents; i++)
          setDofVector(getVecSolInterior(), vec.getDOFVector(i), i, true);

        vecSolAssembly();
      }


      MatNullSpace matNullspace;
      Vec nullspaceBasis;
      if (nullspace.size() > 0 ||
          hasConstantNullspace ||
          constNullspaceComponent.size() > 0)
      {
        TEST_EXIT_DBG(nullspace.size() <= 1)("Not yet implemented!\n");

        if (constNullspaceComponent.size() > 0)
        {
          nullspace.clear();
          SystemVector* basisVec = new SystemVector(vec);
          basisVec->set(0.0);
          for (unsigned int i = 0; i < constNullspaceComponent.size(); i++)
            basisVec->getDOFVector(constNullspaceComponent[i])->set(1.0);

          nullspace.push_back(basisVec);
        }

        if (nullspace.size() > 0)
        {
          VecDuplicate(getVecSolInterior(), &nullspaceBasis);
          setDofVector(nullspaceBasis, *(nullspace[0]), true);

          VecAssemblyBegin(nullspaceBasis);
          VecAssemblyEnd(nullspaceBasis);

          VecNormalize(nullspaceBasis, PETSC_NULL);

          MatNullSpaceCreate(domainComm, (hasConstantNullspace ? PETSC_TRUE : PETSC_FALSE),
                             1, &nullspaceBasis, &matNullspace);

          MatMult(getMatInterior(), nullspaceBasis, getVecSolInterior());
          PetscReal n;
          VecNorm(getVecSolInterior(), NORM_2, &n);
          MSG("NORM IS: %e\n", n);
        }
        else
        {
          MatNullSpaceCreate(domainComm, PETSC_TRUE, 0, PETSC_NULL, &matNullspace);
        }

        MSG("NULLSPACE IS NOT REMOVED!\n");
        // MatSetNullSpace(getMatInterior(), matNullspace);
        // KSPSetNullSpace(kspInterior, matNullspace);

        // === Remove null space, if requested. ===

        if (removeRhsNullspace)
        {
          TEST_EXIT_DBG(coarseSpaceMap.empty())("Not supported!\n");

          MSG("Remove nullspace from rhs vector.\n");
#if (PETSC_VERSION_MINOR >= 5)
          MatNullSpaceRemove(matNullspace, getVecRhsInterior());
#else
          MatNullSpaceRemove(matNullspace, getVecRhsInterior(), PETSC_NULL);
#endif
        }
      }
      else
      {
        TEST_EXIT(removeRhsNullspace == false)
        ("No nullspace provided that should be removed from rhs!\n");
      }

      // PETSc.
      solve(getVecRhsInterior(), getVecSolInterior());


      if (nullspace.size() > 0)
      {
        MatNullSpaceDestroy(&matNullspace);
        VecDestroy(&nullspaceBasis);
      }


      // === Transfere values from PETSc's solution vectors to the DOF vectors. ===
      PetscScalar* vecPointer;
      VecGetArray(getVecSolInterior(), &vecPointer);

      int c = 0;
      for (int component = 0; component < nComponents; component++)
      {
        DOFVector<double>& dv = *(vec.getDOFVector(component));

        DofMap& d = (*interiorMap)[component].getMap();
        for (DofMap::iterator it = d.begin(); it != d.end(); ++it)
          if (it->second.local != -1)
            dv[it->first] = vecPointer[c++];
      }

      VecRestoreArray(getVecSolInterior(), &vecPointer);

      // === Synchronize DOFs at common DOFs, i.e., DOFs that correspond to ===
      // === more than one partition.                                       ===
      meshDistributor->synchVector(vec);
    }


    void PetscSolverGlobalMatrix::solveGlobal(Vec& rhs, Vec& sol)
    {
      Vec tmp;
      if (domainComm.Get_size() == 1)
        createLocalVec(*interiorMap, tmp);
      else
        createVec(*interiorMap, tmp);

      PetscScalar* tmpValues, *rhsValues;
      VecGetArray(tmp, &tmpValues);
      VecGetArray(rhs, &rhsValues);

      for (int i = 0; i < interiorMap->getRankDofs(); i++)
        tmpValues[i] = rhsValues[i];

      VecRestoreArray(rhs, &rhsValues);
      VecRestoreArray(tmp, &tmpValues);

      KSPSolve(kspInterior, tmp, tmp);

      VecGetArray(tmp, &tmpValues);
      VecGetArray(sol, &rhsValues);

      for (int i = 0; i < interiorMap->getRankDofs(); i++)
        rhsValues[i] = tmpValues[i];

      VecRestoreArray(sol, &rhsValues);
      VecRestoreArray(tmp, &tmpValues);

      VecDestroy(&tmp);
    }


    void PetscSolverGlobalMatrix::destroyMatrixData()
    {
      matDestroy();

      exitPreconditioner(pcInterior);
      exitSolver(kspInterior);
    }


    void PetscSolverGlobalMatrix::destroyVectorData()
    {
      vecDestroy();
    }


    void PetscSolverGlobalMatrix::createFieldSplit(PC pc)
    {
      FUNCNAME("PetscSolverGlobalMatrix::createFieldSplit()");

      vector<string> isNames;
      Parameters::get("parallel->solver->is blocks", isNames);

      int nBlocks = isNames.size();
      if (nBlocks == 0)
        return;

      for (int i = 0; i < nBlocks; i++)
      {
        MSG("Create for block %s\n", isNames[i].c_str());

        vector<int> blockComponents;
        Parameters::get("parallel->solver->is block " + std::to_string(i),
                        blockComponents);
        int nComponents = static_cast<int>(blockComponents.size());

        TEST_EXIT(nComponents > 0)("No IS block for block %d defined!\n", i);

        // Check if blocks are continous
        for (int j = 0; j < nComponents; j++)
        {
          TEST_EXIT(blockComponents[j] == blockComponents[0] + j)
          ("Does not yet support not continous IS blocks! Block %s\n",
           isNames[i].c_str());
        }

        createFieldSplit(pc, isNames[i].c_str(), blockComponents);
      }
    }


    void PetscSolverGlobalMatrix::createFieldSplit(PC pc,
        const char* splitName,
        vector<int>& components)
    {
      IS is;
      interiorMap->createIndexSet(is, components[0], components.size());
      PCFieldSplitSetIS(pc, splitName, is);
      ISDestroy(&is);
    }


    void PetscSolverGlobalMatrix::extractVectorComponent(Vec input, int i, Vec* output, int numberOfComponents)
    {
      IS is;
      interiorMap->createIndexSet(is, i, numberOfComponents);
      VecGetSubVector(input, is, output);
      ISDestroy(&is);
    }

    void PetscSolverGlobalMatrix::extractMatrixComponent(Mat input, int startRow, int numberOfRows, int startCol, int numberOfCols, Mat* output)
    {
      IS isrow, iscol;
      interiorMap->createIndexSet(isrow, startRow, numberOfRows);
      interiorMap->createIndexSet(iscol, startCol, numberOfCols);
      MatGetSubMatrix(input, isrow, iscol, MAT_INITIAL_MATRIX, output);
      ISDestroy(&iscol);
      ISDestroy(&isrow);
    }


    void PetscSolverGlobalMatrix::initSolver(KSP& ksp)
    {
      KSPCreate(domainComm, &ksp);
#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(ksp, getMatInterior(), getMatInterior());
#else
      KSPSetOperators(ksp, getMatInterior(), getMatInterior(), SAME_NONZERO_PATTERN);
#endif
      KSPSetTolerances(ksp, 0.0, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT);
      KSPSetType(ksp, KSPBCGS);
      KSPSetOptionsPrefix(ksp, kspPrefix.c_str());
      MatSetOptionsPrefix(getMatInterior(), kspPrefix.c_str());
      KSPSetFromOptions(ksp);

      // Do not delete the solution vector, use it for the initial guess.
      if (!zeroStartVector)
        KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    }


    void PetscSolverGlobalMatrix::exitSolver(KSP& ksp)
    {
      KSPDestroy(&ksp);
    }


    void PetscSolverGlobalMatrix::initPreconditioner(PC pc)
    {
      PCSetFromOptions(pc);
      createFieldSplit(pc);
    }


    void PetscSolverGlobalMatrix::exitPreconditioner(PC pc)
    { }


    /// called by \ref fillPetscMatrix
    void PetscSolverGlobalMatrix::setDofMatrix(DOFMatrix* seqMat,
        int rowComp, int colComp)
    {
      FUNCNAME("PetscSolverGlobalMatrix::setDofMatrix()");

      TEST_EXIT(seqMat)("No DOFMatrix!\n");

      using mtl::tag::row;
      using mtl::tag::nz;
      using mtl::begin;
      using mtl::end;
      namespace traits = mtl::traits;
      typedef DOFMatrix::base_matrix_type Matrix;

      traits::col<Matrix>::type col(seqMat->getBaseMatrix());
      traits::const_value<Matrix>::type value(seqMat->getBaseMatrix());

      typedef traits::range_generator<row, Matrix>::type cursor_type;
      typedef traits::range_generator<nz, cursor_type>::type icursor_type;

      vector<int> cols;
      vector<double> values;
      cols.reserve(300);
      values.reserve(300);

      vector<int> globalCols;

      // Get periodic mapping object
      PeriodicMap& perMap = meshDistributor->getPeriodicMap();
      std::set<DegreeOfFreedom>& dirichletRows = seqMat->getDirichletRows();

      const FiniteElemSpace* rowFe = seqMat->getRowFeSpace();
      const FiniteElemSpace* colFe = seqMat->getColFeSpace();

      // === Traverse all rows of the DOF matrix and insert row wise the values ===
      // === to the PETSc matrix.                                               ===

      for (cursor_type cursor = begin<row>(seqMat->getBaseMatrix()),
           cend = end<row>(seqMat->getBaseMatrix()); cursor != cend; ++cursor)
      {
        // Global index of the current row DOF.
        MultiIndex rowMultiIndex;
        if ((*interiorMap)[rowComp].find(cursor.value(), rowMultiIndex) == false)
          continue;

        int globalRowDof = rowMultiIndex.global;

        // Test if the current row DOF is a periodic DOF.
        bool periodicRow = perMap.isPeriodic(rowFe, globalRowDof);

        // Dirichlet rows can be set only be the owner ranks.
        if (dirichletRows.count(cursor.value()) && !((*interiorMap)[rowComp].isRankDof(cursor.value())))
          continue;

        if (!periodicRow)
        {
          // === Row DOF index is not periodic. ===

          // Get PETSc's mat row index.
          int rowIndex = interiorMap->getMatIndex(rowComp, globalRowDof);

          cols.clear();
          values.clear();

          for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor);
               icursor != icend; ++icursor)
          {

            // Global index of the current column index.
            MultiIndex colMultiIndex;
            if ((*interiorMap)[colComp].find(col(*icursor), colMultiIndex) == false)
              continue;

            int globalColDof = colMultiIndex.global;
            // Test if the current col dof is a periodic dof.
            bool periodicCol = perMap.isPeriodic(colFe, globalColDof);
            // Get PETSc's mat col index.
            int colIndex = interiorMap->getMatIndex(colComp, globalColDof);

            // Ignore all zero entries, expect it is a diagonal entry.
            if (value(*icursor) == 0.0 && rowIndex != colIndex)
              continue;

            if (!periodicCol)
            {
              // Calculate the exact position of the column index in the PETSc matrix.
              cols.push_back(colIndex);
              values.push_back(value(*icursor));
            }
            else
            {
              // === Row index is not periodic, but column index is. ===

              // Create set of all periodic associations of the column index.
              std::set<int> perAsc;
              perMap.fillAssociations(colFe, globalColDof,
                                      meshDistributor->getElementObjectDb(), perAsc);

              // Scale value to the number of periodic associations of the column index.
              double scaledValue =
                value(*icursor) * pow(0.5, static_cast<double>(perAsc.size()));


              // === Create set of all matrix column indices due to the periodic ===
              // === associations of the column DOF index.                       ===

              vector<int> newCols;
              perMap.mapDof(colFe, globalColDof, perAsc, newCols);
              for (unsigned int i = 0; i < newCols.size(); i++)
              {
                cols.push_back(interiorMap->getMatIndex(colComp, newCols[i]));
                values.push_back(scaledValue);
              }
            }
          }

          MatSetValues(getMatInterior(), 1, &rowIndex, cols.size(),
                       &(cols[0]), &(values[0]), ADD_VALUES);
        }
        else
        {
          // === Row DOF index is periodic. ===

          // Because this row is periodic, we will have to add the entries of this
          // matrix row to multiple rows. The following maps store to each row an
          // array of column indices and values of the entries that must be added to
          // the PETSc matrix.
          map<int, vector<int>> colsMap;
          map<int, vector<double>> valsMap;

          // Traverse all column entries.
          for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor);
               icursor != icend; ++icursor)
          {
            // Global index of the current column index.
            int globalColDof = (*interiorMap)[colComp][col(*icursor)].global;

            // Ignore all zero entries, expect it is a diagonal entry.
            if (value(*icursor) == 0.0 && globalRowDof != globalColDof)
              continue;

            // === Add all periodic associations of both, the row and the column ===
            // === indices to the set perAsc.                                    ===

            std::set<int> perAsc;
            perMap.fillAssociations(colFe, globalColDof,
                                    meshDistributor->getElementObjectDb(), perAsc);
            perMap.fillAssociations(rowFe, globalRowDof,
                                    meshDistributor->getElementObjectDb(), perAsc);

            // Scale the value with respect to the number of periodic associations.
            double scaledValue =
              value(*icursor) * pow(0.5, static_cast<double>(perAsc.size()));


            // === Create all matrix entries with respect to the periodic  ===
            // === associations of the row and column indices.             ===

            vector<pair<int, int>> entry;
            perMap.mapDof(rowFe, colFe, make_pair(globalRowDof, globalColDof),
                          perAsc, entry);

            // === Translate the matrix entries to PETSc's matrix.

            for (unsigned int i = 0; i < entry.size(); i++)
            {
              int rowIdx = interiorMap->getMatIndex(rowComp, entry[i].first);
              int colIdx = interiorMap->getMatIndex(colComp, entry[i].second);

              colsMap[rowIdx].push_back(colIdx);
              valsMap[rowIdx].push_back(scaledValue);
            }
          }


          // === Finally, add all periodic rows to the PETSc matrix. ===

          for (map<int, vector<int>>::iterator rowIt = colsMap.begin();
               rowIt != colsMap.end(); ++rowIt)
          {
            TEST_EXIT_DBG(rowIt->second.size() == valsMap[rowIt->first].size())
            ("Should not happen!\n");

            int rowIndex = rowIt->first;

            MatSetValues(getMatInterior(), 1, &rowIndex, rowIt->second.size(),
                         &(rowIt->second[0]), &(valsMap[rowIt->first][0]), ADD_VALUES);
          }
        }
      }
    }


    /// called by \ref fillPetscRhs
    void PetscSolverGlobalMatrix::setDofVector(Vec vecInterior,
        Vec vecCoarse,
        DOFVector<double>* vec,
        int rowComp,
        bool rankOnly)
    {
      FUNCNAME_DBG("PetscSolverGlobalMatrix::setDofVector()");

      const FiniteElemSpace* feSpace = vec->getFeSpace();
      PeriodicMap& perMap = meshDistributor->getPeriodicMap();

      ParallelDofMapping* rowCoarseSpace =
        (coarseSpaceMap.size() ? coarseSpaceMap[rowComp] : NULL);

      map<DegreeOfFreedom, double>& dirichletValues = vec->getDirichletValues();

      // Traverse all used DOFs in the dof vector.
      DOFVector<double>::Iterator dofIt(vec, USED_DOFS);
      for (dofIt.reset(); !dofIt.end(); ++dofIt)
      {

        DegreeOfFreedom dof = dofIt.getDOFIndex();

        if (rankOnly && !(*interiorMap)[rowComp].isRankDof(dof))
          continue;

        bool isCoarseDof = isCoarseSpace(rowComp, dof);

        // Dirichlet rows can be set only be the owner ranks.
        if (dirichletValues.count(dof))
        {
          if ((!isCoarseDof && !((*interiorMap)[rowComp].isRankDof(dof))) ||
              (isCoarseDof && !((*rowCoarseSpace)[rowComp].isRankDof(dof))))
            continue;
        }

        if (isCoarseDof)
        {
          TEST_EXIT_DBG(vecCoarse != PETSC_NULL)("vecCoarse not set! Should not happen!\n");

          int index = rowCoarseSpace->getMatIndex(rowComp, dof);
          VecSetValue(vecCoarse, index, *dofIt, ADD_VALUES);
        }
        else
        {
          if ((*interiorMap)[rowComp].isSet(dof) == false)
            continue;

          // Calculate global row index of the DOF.
          DegreeOfFreedom globalRowDof = (*interiorMap)[rowComp][dof].global;

          // Get PETSc's mat index of the row DOF.
          int index = 0;
          if (interiorMap->isMatIndexFromGlobal())
            index = interiorMap->getMatIndex(rowComp, globalRowDof) + rStartInterior;
          else
            index =
              interiorMap->getMatIndex(rowComp, dof) + rStartInterior;

          if (perMap.isPeriodic(feSpace, globalRowDof))
          {
            std::set<int>& perAsc = perMap.getAssociations(feSpace, globalRowDof);
            double value = *dofIt / (perAsc.size() + 1.0);
            VecSetValue(vecInterior, index, value, ADD_VALUES);

            for (std::set<int>::iterator perIt = perAsc.begin();
                 perIt != perAsc.end(); ++perIt)
            {
              int mappedDof = perMap.map(feSpace, *perIt, globalRowDof);
              int mappedIndex = interiorMap->getMatIndex(rowComp, mappedDof);

              VecSetValue(vecInterior, mappedIndex, value, ADD_VALUES);
            }
          }
          else
          {
            // The DOF index is not periodic.
            VecSetValue(vecInterior, index, *dofIt, ADD_VALUES);
          }
        }
      }
    }


    PetscSolver* PetscSolverGlobalMatrix::createSubSolver(int component,
        string kspPrefix)
    {
      vector<const FiniteElemSpace*> fe;
      fe.push_back(componentSpaces[component]);

      PetscSolver* subSolver = new PetscSolverGlobalMatrix("");
      subSolver->setKspPrefix(kspPrefix);
      subSolver->setMeshDistributor(meshDistributor, 0);
      subSolver->init(fe, fe);

      ParallelDofMapping& subDofMap = subSolver->getDofMap();
      subDofMap[0] = dofMap[component];
      subDofMap.update();

      return subSolver;
    }


    void PetscSolverGlobalMatrix::setConstantNullSpace(KSP ksp,
        int constFeSpace,
        bool test)
    {
      FUNCNAME("PetscSolverGlobalMatrix::setConstantNullSpace()");

      Vec nullSpaceBasis;
      VecDuplicate(getVecSolInterior(), &nullSpaceBasis);

      SystemVector basisVec("tmp", componentSpaces, componentSpaces.size(), true);
      basisVec.set(0.0);
      basisVec.getDOFVector(constFeSpace)->set(1.0);

      setDofVector(nullSpaceBasis, basisVec, true);
      VecAssemblyBegin(nullSpaceBasis);
      VecAssemblyEnd(nullSpaceBasis);
      VecNormalize(nullSpaceBasis, PETSC_NULL);

      if (test)
      {
        Vec tmp;
        MatGetVecs(getMatInterior(), &tmp, PETSC_NULL);
        MatMult(getMatInterior(), nullSpaceBasis, tmp);
        PetscReal n;
        VecNorm(tmp, NORM_2, &n);
        MSG("NORM IS: %e\n", n);
        VecDestroy(&tmp);
      }


      MatNullSpace matNullSpace;
      MatNullSpaceCreate(domainComm, PETSC_FALSE, 1, &nullSpaceBasis, &matNullSpace);
      KSPSetNullSpace(ksp, matNullSpace);
      MatNullSpaceDestroy(&matNullSpace);

      VecDestroy(&nullSpaceBasis);
    }


    void PetscSolverGlobalMatrix::setConstantNullSpace(KSP ksp)
    {
      MatNullSpace matNullSpace;
      MatNullSpaceCreate(domainComm, PETSC_TRUE, 0, PETSC_NULL, &matNullSpace);
      KSPSetNullSpace(ksp, matNullSpace);
      MatNullSpaceDestroy(&matNullSpace);
    }

  }
} // end namespace Parallel, AMDiS
