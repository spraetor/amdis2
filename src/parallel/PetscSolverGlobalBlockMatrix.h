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



/** \file PetscSolverGlobalBlockMatrix.h */

#ifndef AMDIS_PETSC_SOLVER_GLOBAL_BLOCK_MATRIX_H
#define AMDIS_PETSC_SOLVER_GLOBAL_BLOCK_MATRIX_H

#include "AMDiS_fwd.h"
#include "parallel/PetscSolver.h"

namespace AMDiS
{
  namespace Parallel
  {

    /** \ingroup Solver
     *
     * \brief
     * PETSc solver which creates a globally distributed (nested) matrix.
     */
    class PetscSolverGlobalBlockMatrix : public PetscSolver
    {
    public:
      /// Creator class
      class Creator : public LinearSolverCreator
      {
      public:
        virtual ~Creator() {}

        /// Returns a new PetscSolver object.
        LinearSolverInterface* create()
        {
          return new PetscSolverGlobalBlockMatrix(this->name);
        }
      };

      PetscSolverGlobalBlockMatrix(std::string name)
        : PetscSolver(name),
          nComponents(0),
          nBlocks(-1)
      {}

      void fillPetscMatrix(Matrix<DOFMatrix*>* mat);

      void fillPetscRhs(SystemVector* vec);

      void solvePetscMatrix(SystemVector& vec, AdaptInfo* adaptInfo);

      void destroyMatrixData();

      void destroyVectorData();

    protected:
      /// Takes a DOF matrix and sends the values to the global PETSc matrix.
      void setDofMatrix(Mat& petscMat, DOFMatrix* mat,
                        int dispRowBlock, int dispColBlock);

      /// Takes a DOF vector and sends its values to a given PETSc vector.
      void setDofVector(Vec& petscVec, DOFVector<double>* vec);

      virtual void initSolver(KSP& ksp);

      virtual void exitSolver(KSP ksp);

      virtual void initPreconditioner(PC pc);

      virtual void exitPreconditioner(PC pc);

    protected:
      std::vector<Mat> nestMat;

      std::vector<Vec> nestVec;

      Vec petscSolVec;

      /// Number of components (= number of unknowns in the PDE)
      int nComponents;

      /// Number of blocks for the solver, must be 1 <= nBlocks <= nComponents
      int nBlocks;

      /// Maps to each component number the block number the component is in.
      std::map<int, int> componentInBlock;
    };

  } // end namespace Parallel

} // end namespace AMDiS

#endif

