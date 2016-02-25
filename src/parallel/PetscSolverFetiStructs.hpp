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



/** \file PetscSolverFetiStructs.h */

#ifndef AMDIS_PETSC_SOLVER_FETI_STRUCTS_H
#define AMDIS_PETSC_SOLVER_FETI_STRUCTS_H

#include <map>
#include "parallel/PetscSolver.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    class PetscSolverFeti;


    enum FetiSolverType
    {
      // Standard exakt FETI-DP system
      EXACT,
      // Inexact FETI-DP
      INEXACT,
      // Inexact reduced FETI-DP
      INEXACT_REDUCED
    };

    /** \brief
    * This structure is used when defining the MatShell operation for solving
    * primal schur complement. \ref petscMultMatSchurPrimal
    */
    struct SchurPrimalData
    {
      /// Temporal vector on the B variables.
      Vec tmp_vec_b;

      /// Temporal vecor in the primal variables.
      Vec tmp_vec_primal;

      PetscSolver* subSolver;
    };


    /** \brief
    *
    */
    struct SchurPrimalAugmentedData
    {
      /// Temporal vectors on the B variables.
      Vec tmp_vec_b0, tmp_vec_b1;

      Vec tmp_vec_primal;

      Vec tmp_vec_lagrange;

      Mat* mat_lagrange;

      Mat* mat_augmented_lagrange;

      PetscSolver* subSolver;

      bool nestedVec;
    };


    /** \brief
    * This structure is used when defining the FETI-DP operator for solving
    * the system matrix reduced to the Lagrange multipliers.
    * \ref petscMultMatFeti
    */
    struct FetiData
    {
      /// Matrix of Lagrange variables.
      Mat* mat_lagrange;

      ///
      Mat* mat_augmented_lagrange;

      /// Temporal vectors on the B variables.
      Vec tmp_vec_b0, tmp_vec_b1;

      /// Temporal vector on the primal variables.
      Vec tmp_vec_primal0, tmp_vec_primal1;

      /// Temporal vector on the lagrange variables.
      Vec tmp_vec_lagrange;

      Vec tmp_vec_interface;

      PetscSolver* subSolver;

      /// Pointer to the solver of the schur complement on the primal variables.
      KSP* ksp_schur_primal;
    };


    struct FetiInexactData
    {
      Mat* matBB, *matBPi, *matPiB, *matPiPi;

      Mat* mat_lagrange;

      Vec tmp_vec_b0, tmp_vec_b1;
    };


    struct FetiInexactPreconData
    {
      KSP ksp_schur;

      KSP ksp_interior;

      KSP ksp_pc_feti;

      PC pc_feti;

      Mat* matPiB, *matBPi;

      Mat* mat_lagrange;

      Vec tmp_vec_b0;
    };


    struct FetiDirichletPreconData
    {
      /// Matrix of scaled Lagrange variables.
      Mat* mat_lagrange_scaled;

      Mat* mat_interior_interior, *mat_duals_duals;

      Mat* mat_interior_duals, *mat_duals_interior;

      /// Pointer to the solver for \ref PetscSolverFeti::mat_bb.
      KSP* ksp_interior;

      /// Temporal vector on the B variables.
      Vec tmp_vec_b;

      /// Temporal vector on the dual variables.
      Vec tmp_vec_duals0, tmp_vec_duals1;

      /// Temporal vector on the interior variables.
      Vec tmp_vec_interior;

      std::map<int, int> localToDualMap;
    };


    struct FetiLumpedPreconData
    {
      /// Matrix of scaled Lagrange variables.
      Mat* mat_lagrange_scaled;

      Mat* mat_duals_duals;

      /// Temporal vector on the B variables.
      Vec tmp_vec_b0;

      /// Temporal vector on the dual variables.
      Vec tmp_vec_duals0, tmp_vec_duals1;

      std::map<int, int> localToDualMap;
    };


    struct FetiInterfaceLumpedPreconData : public FetiLumpedPreconData
    {
      /// Temporal vectors on the B variables.
      Vec tmp_vec_b1;

      PetscSolver* subSolver;

      Vec tmp_primal;

      KSP ksp_mass;
    };


    struct FetiKspData
    {
      Vec draft;
    };

    typedef enum
    {
      FETI_NONE = 0,
      FETI_DIRICHLET = 1,
      FETI_LUMPED = 2
    } FetiPreconditioner;
  }
}

#endif
