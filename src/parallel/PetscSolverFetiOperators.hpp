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



/** \file PetscSolverFeti.h */

#ifndef AMDIS_PETSC_SOLVER_FETI_OPERATORS_H
#define AMDIS_PETSC_SOLVER_FETI_OPERATORS_H

#include <mpi.h>
#include <petsc.h>

namespace AMDiS
{
  namespace Parallel
  {
    void copyGlobalLocal(Vec globalB, Vec localB);

    ///
    int petscMultMatSchurPrimal(Mat mat, Vec x, Vec y);

    ///
    int petscMultMatSchurPrimalAugmented(Mat mat, Vec x, Vec y);

    /// FETI-DP operator
    int petscMultMatFeti(Mat mat, Vec x, Vec y);

    /// Inexact FETI-DP operator
    int petscMultMatFetiInexact(Mat mat, Vec x, Vec y);

    /// Inexact FETI-DP preconditioner
    PetscErrorCode pcInexactFetiShell(PC pc, Vec x, Vec y);

    /// FETI-DP operator with augmented Lagrange constraints
    int petscMultMatFetiAugmented(Mat mat, Vec x, Vec y);

    /// FETI-DP operator used for Stokes like problems.
    int petscMultMatFetiInterface(Mat mat, Vec x, Vec y);

    /// y = PC * x
    PetscErrorCode petscApplyFetiDirichletPrecon(PC pc, Vec x, Vec y);

    /// y = PC * x
    PetscErrorCode petscApplyFetiLumpedPrecon(PC pc, Vec x, Vec y);

    PetscErrorCode petscApplyFetiInterfaceLumpedPrecon(PC pc, Vec x, Vec y);
  }
}

#endif
