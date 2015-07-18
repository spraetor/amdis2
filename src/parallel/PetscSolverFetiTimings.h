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



/** \file PetscSolverFetiTimings.h */


#ifndef AMDIS_PETSC_SOLVER_FETI_TIMINGS_H
#define AMDIS_PETSC_SOLVER_FETI_TIMINGS_H

namespace AMDiS
{
  namespace Parallel
  {

    class FetiTimings {
    private:
      FetiTimings() {}

    public:
      static void reset();

    public:
      static double fetiSolve;
      static double fetiSolve01;
      static double fetiSolve02;

      static double fetiPreconditioner;
    };
    
  } // end namespace Parallel
} // end namespace AMDiS

#endif
