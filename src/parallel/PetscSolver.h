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



/** \file PetscSolver.h */

#ifndef AMDIS_PETSC_SOLVER_H
#define AMDIS_PETSC_SOLVER_H

#include <mpi.h>
#include <petsc.h>
#include <petscsys.h>
#include <petscao.h>
#include <petscksp.h>

#include "AMDiS_fwd.h"
#include "parallel/ParallelCoarseSpaceSolver.h"

namespace AMDiS
{
  namespace Parallel
  {

    /** \ingroup Solver
     *
     * \brief
    * Create an abstract interface to an arbitrary PETSc solver. This class is
    * based on \ref ParallelCoarseSpaceSolver to support for solvers which make
    * use of a coarse grid problem.
    */
    class PetscSolver : public ParallelCoarseSpaceSolver
    {
    public: // typedefs
      typedef ParallelCoarseSpaceSolver super;

    public: // methods
      PetscSolver(std::string name);

      virtual ~PetscSolver();

      /// calls super::init() and sets boundaryDofRequirement for meshdistributor
      virtual void init(std::vector<const FiniteElemSpace*>& componentSpaces,
                        std::vector<const FiniteElemSpace*>& feSpaces,
                        bool createGlobalMapping = true);

      /** \brief
      * Create a PETSc matrix. The given DOF matrices are used to create the nnz
      * structure of the PETSc matrix and the values are transfered to it.
      *
      * \param[in] mat
      */
      virtual void fillPetscMatrix(Matrix<DOFMatrix*>* mat) = 0;

      /// This function is just a small wrapper that creates a 1x1 matrix that
      /// contains exactly one DOFMatrix and than calls \ref fillPetscMatrix
      void fillPetscMatrix(DOFMatrix* mat);

      /** \brief
      * Create a PETSc vector and fills it with the rhs values of the system.
      *
      * \param[in] vec
      */
      virtual void fillPetscRhs(SystemVector* vec) = 0;

      /// Use PETSc to solve the linear system of equations
      virtual void solvePetscMatrix(SystemVector& vec, AdaptInfo* adaptInfo) = 0;

      virtual void solve(Vec& rhs, Vec& sol);

      virtual void solveGlobal(Vec& rhs, Vec& sol);

      /// Destroys all matrix data structures.
      virtual void destroyMatrixData() = 0;

      /// Detroys all vector data structures.
      virtual void destroyVectorData() = 0;

      virtual Flag getBoundaryDofRequirement()
      {
        return 0;
      }

      KSP getSolver()
      {
        return kspInterior;
      }

      PC getPc()
      {
        return pcInterior;
      }

      void setKspPrefix(std::string s)
      {
        kspPrefix = s;
      }

      void setRemoveRhsNullspace(bool b)
      {
        removeRhsNullspace = b;
      }

      void setSymmetric(bool b)
      {
        isSymmetric = b;
      }

      /// Adds a new vector to the basis of the operator's nullspace.
      void addNullspaceVector(SystemVector* vec)
      {
        nullspace.push_back(vec);
      }

      /// Sets the nullspace to be constant for some specific components.
      void setConstantNullspace(std::vector<int>& components)
      {
        constNullspaceComponent = components;
      }

      /// Sets the nullspace to be constant for a specific component.
      void setConstantNullspace(int component)
      {
        constNullspaceComponent.clear();
        constNullspaceComponent.push_back(component);
      }

      /// Informs the solver whether is has to handle dirichlet rows or not.
      void setHandleDirichletRows(bool b)
      {
        handleDirichletRows = b;
      }

    protected:
      /// Implementation of \ref LinearSolverInterface::solveLinearSystem()
      int solveLinearSystem(const SolverMatrix<Matrix<DOFMatrix*>>& A,
                            SystemVector& x,
                            SystemVector& b,
                            bool createMatrixData,
                            bool storeMatrixData);

      /** \brief
      * Copies between to PETSc vectors by using different index sets for the
      * origin and the destination vectors.
      *
      * \param[in]   originVec    The PETSc vector from which we copy from.
      * \param[out]  destVec      The PETSc vector we copy too.
      * \param[in]   originIndex  Set of global indices referring to the
      *                           origin vector.
      * \param[in]   destIndex    Set of global indices referring to the
      *                           destination vector.
      */
      void copyVec(Vec& originVec, Vec& destVec,
                   std::vector<int>& originIndex, std::vector<int>& destIndex);

      /// Run test, if matrix is symmetric.
      bool testMatrixSymmetric(Mat mat, bool advancedTest = false);

    protected:
      /// PETSc solver object
      KSP kspInterior;

      /// PETSc preconditioner object
      PC pcInterior;

      /// A set of vectors that span the null space of the operator.
      std::vector<SystemVector*> nullspace;

      /// KSP database prefix
      std::string kspPrefix;

      /// If true, the constant null space is projected out of the RHS vector. It
      /// depends on the specific PETSc solver if it considers this value.
      bool removeRhsNullspace;

      bool hasConstantNullspace;

      bool isSymmetric;

      /// If true, dirichlet rows are handled by the solver correspondently. To
      /// set this value to false makes only sense, of this solver is just used
      /// as a subsolver and the main solver above alread handles dirichlet rows
      /// in some way.
      bool handleDirichletRows;

      std::vector<int> constNullspaceComponent;
    };
  } // end namespace Parallel
} // end namespace AMDiS

#endif
