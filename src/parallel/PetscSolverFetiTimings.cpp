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


#include "parallel/PetscSolverFetiTimings.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    double FetiTimings::fetiSolve = 0.0;
    double FetiTimings::fetiSolve01 = 0.0;
    double FetiTimings::fetiSolve02 = 0.0;
    double FetiTimings::fetiPreconditioner = 0.0;


    void FetiTimings::reset()
    {
      fetiSolve = 0.0;
      fetiSolve01 = 0.0;
      fetiSolve02 = 0.0;
      fetiPreconditioner = 0.0;
    }

  }
}
