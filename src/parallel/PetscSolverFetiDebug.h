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



/** \file PetscSolverFetiDebug.h */

#ifndef AMDIS_PETSC_SOLVER_FETI_DEBUG_H
#define AMDIS_PETSC_SOLVER_FETI_DEBUG_H

namespace AMDiS
{
  namespace Parallel 
  {

    /** \brief
    * This class collects several static functions which may be used to debug
    * FETI-DP solver. All functions should be used only in debug mode and
    * never in productive runs as they may be arbitrarly slow.
    *
    * This class must be a friend class of \ref PetscSolverFeti as it must have
    * access to all internal data.
    */
    class PetscSolverFetiDebug
    {
    public:
      /** \brief
      * Is used to test the null space in stokes mode. The function reads the
      * null space basis from FETI mat object and test by simple matrix 
      * multiplication, whether the vector is a member of the null space.
      * Furthermore, the function creates the null space on the non reduced
      * form of the FETI-DP system and make there the same test. The resulting
      * test vector is also written to "nullspace.vtu" file.
      *
      *
      * \param[in]  feti    FETI-DP solver
      * \param[in]  vec     Is just used as a temporary template vector.
      */
      static void debugNullSpace(PetscSolverFeti &feti, SystemVector &vec);

      /** \brief
      * Creates a parallel distributed matrix of the local (B-index) matrices. 
      * Note that the matrix is still a block matrix and there are no off 
      * diagonal values.
      *
      * \param[in]  feti     FETI-DP solver
      * \param[out] mat      The parallel distributed matrix will be created
      *                      on this variable.
      */
      static void createInteriorMat(PetscSolverFeti &feti, Mat &mat);

      /** \brief
      * Creates the non reduced FETI-DP system (also in stokes mode) as a 
      * nested matrix.
      *
      * \param[in]  feti     FETI-DP solver
      * \param[out] mat      The resulting nested matrix.
      */
      static void createNestedFetiMat(PetscSolverFeti &feti, Mat &mat);

      /** \brief
      * Creates the reduced FETI-DP operator in an explicitly given matrix
      * form (also in stokes mode). The resulting matrix is dense.
      *
      * \param[in]  feti         FETI-DP solver
      * \param[in]  fetiMat      Matrix object representing the implicit
      *                          FETI-DP operator
      * \param[out] explicitMat  Explicit matrix formulation
      * \param[out] nnzCounter   Stores the number of non zero elements in the
      *                          explicit matrix formulation
      */
      static void createExplicitFetiMat(PetscSolverFeti &feti, 
					Mat fetiMat, 
					Mat &explicitmat, 
					int &nnzCounter);

      /** \brief
      * Create a non nested parallel distributed vector from a nested vector of
      * two parallel distributed vectors. Is used to change the structure of
      * solution vectors in stokes mode, which always have two components.
      *
      * \param[in]  feti         FETI-DP solver
      * \param[in]  nestedVec    Nested PETSc vector with two sub vectors.
      * \param[out] explicitVec  The output vector.
      */
      static void createExplicitVec(PetscSolverFeti &feti, 
				    Vec nestedVec, 
				    Vec &explicitVec);

      /** \brief
      * Checks for the init file parameter "parallel->debug->write fety system"
      * and if it is true, writes the null space basis vector to the PETSc binary
      * file "nullspace.vec".
      *
      * \param[in] feti            FETI-DP solver
      * \param[in] nullSpaceBasis  Vector representing the basis of the null
      *                            space of the FETI-DP system.
      */
      static void writeNullSpace(PetscSolverFeti &feti, 
				Vec nullSpaceBasis);

      /** \brief
      * Checks for the init file parameter "parallel->debug->write fety system"
      * and if it is true, some debug test of the FETI-DP system are made:
      *   - symmetry test of the FETI-DP operator (by assembling it in an
      *     explicit formulation)
      *   - write explicit FETI-DP operator to "feti.mat"
      *   - write rhs to "feti_rhs.vec"
      *
      * \param[in] feti       FETI-DP solver
      * \param[in] dbgRhsVec  Vector representing the rhs vector of the 
      *                       FETI-DP system.
      */
      static void debugFeti(PetscSolverFeti &feti, 
			    Vec dbgRhsVec);

      /** \brief
      * This functions check a PETSc matrix for zero rows. Each zero row index
      * is printed to the screen.
      *
      * \param[in]  mat    PETSc matrix which should be checked for zero rows.
      *
      * \return     int    Number of zero rows in matrix;
      */
      static int testZeroRows(Mat mat);

      /// Write files with information about the primal coarner nodes.
      static void writePrimalFiles(PetscSolverFeti &feti);
    };
  }
}

#endif
