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



/** \file PetscHelper.h */

#ifndef AMDIS_PETSCHELPER_H
#define AMDIS_PETSCHELPER_H

#include <mpi.h>
#include <map>
#include <vector>
#include <petsc.h>
#include "AMDiS_fwd.h"

namespace AMDiS
{
  namespace Parallel
  {

    /** \brief
    * In this namespace, we collect several auxiliary functions for using
    * PETSc in AMDiS. Many of these function may be replaced by new PETSc
    * function in upcoming versions.
    */
    namespace petsc_helper
    {

      /// Defines a PETSc matrix column wise
      typedef std::pair<std::vector<int>, std::vector<double>> SparseCol;
      typedef std::map<int, SparseCol> PetscMatCol;

      /** \brief
      * Returns for a distributed matrix on each rank the local matrix in a
      * sparce column format.
      *
      * \param[in]   mat     PETSc distributerd matrix.
      * \param[out]  matCol  The sparse column represenation of the local matrix.
      */
      void getMatLocalColumn(Mat mat, PetscMatCol& matCol);

      /** \brief
      * Set a local column vector in a distributed matrix.
      *
      * \param[out]  mat     Distributed matrix.
      * \param[in]   column  Column index.
      * \param[in]   vec     Column vector.
      */
      void setMatLocalColumn(Mat mat, int column, Vec vec);

      /** \brief
      * Create a local PETSc vector representing the column of a matrix
      * stored in \ref PetscMatCol type format.
      *
      * \param[in]   column   Sparse column representation.
      * \param[out]  vec      Vector representing one column of the matrix.
      */
      void getColumnVec(const SparseCol& matCol, Vec vec);

      /** \brief
      * Computes the matrix matrix product inv(A) B = C. Matrices B and C
      * are distributed matrices. Matrix A is a local matrix on each MPI
      * task. The overall number of rows of local matrices A must be the
      * number of distriubted rows in B.
      *
      * \param[in]   mpiComm  MPI Communicator object (must fit with ksp)
      * \param[in]   ksp      inv(A) matrix given by a PETSc solver object.
      * \param[in]   mat0     matrix B
      * \param[out]  mat1     resulting matrix C, is created inside the function
      */
      void blockMatMatSolve(MPI::Intracomm mpiComm,
                            KSP ksp,
                            Mat mat0,
                            Mat& mat1);

      /** \brief
      * Converts a 2x2 nested matrix to a MATAIJ matrix (thus not nested).
      *
      * \param[in]  matNest  nested input matrix
      * \param[out] mat      matrix of type MATAIJ, created inside this function.
      */
      void matNestConvert(Mat matNest, Mat& mat);

      void setSolverWithLu(KSP ksp,
                           const char* kspPrefix,
                           KSPType kspType,
                           PCType pcType,
                           const MatSolverPackage matSolverPackage,
                           PetscReal rtol = PETSC_DEFAULT,
                           PetscReal atol = PETSC_DEFAULT,
                           PetscInt maxIt = PETSC_DEFAULT);

      void setSolver(KSP ksp,
                     const char* kspPrefix,
                     KSPType kspType,
                     PCType pcType,
                     PetscReal rtol = PETSC_DEFAULT,
                     PetscReal atol = PETSC_DEFAULT,
                     PetscInt maxIt = PETSC_DEFAULT);

      void createSolver(MPI::Intracomm comm, KSP& ksp, Mat m, std::string kspPrefix = "", int info = 0);

    } // end namespace petsc_helper
  } // end namespace Parallel
} // end namespace AMDiS

#endif
