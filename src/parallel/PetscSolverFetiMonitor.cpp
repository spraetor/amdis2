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


#include "parallel/PetscSolverFetiMonitor.hpp"
#include "parallel/PetscSolverFetiStructs.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    PetscErrorCode KSPMonitorFetiStokes(KSP ksp, PetscInt n, PetscReal rnorm, void* data)
    {
      Vec Br,v,w;
      VecDuplicate(static_cast<FetiKspData*>(data)->draft, &v);
      VecDuplicate(static_cast<FetiKspData*>(data)->draft, &w);

      KSPBuildResidual(ksp, v, w, &Br);

      Vec nest0, nest1;
      VecNestGetSubVec(Br, 0, &nest0);
      VecNestGetSubVec(Br, 1, &nest1);

      PetscScalar norm, norm0, norm1;
      VecNorm(Br, NORM_2, &norm);
      VecNorm(nest0, NORM_2, &norm0);
      VecNorm(nest1, NORM_2, &norm1);

      VecDestroy(&v);
      VecDestroy(&w);

      PetscPrintf(PETSC_COMM_WORLD, "%3D KSP residual norm %1.12e [ %1.12e %1.12e ] and preconditioned norm [%1.12e]\n",
                  n, norm, norm0, norm1, rnorm);

      return 0;
    }

  }
}
