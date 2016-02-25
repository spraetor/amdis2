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

/** \file MtlFetiPrimalSolver.h */

#ifndef AMDIS_MTL_FETI_PRIMAL_SOLVER_H
#define AMDIS_MTL_FETI_PRIMAL_SOLVER_H


#include "AMDiS_fwd.h"
#include "LinearSolverInterface.h"
#include "parallel/PetscSolver.hpp"


namespace AMDiS
{
  namespace Parallel
  {

    /// Solver for the primal nodes based on MTL
    class MtlFetiPrimalSolver
      : public PetscSolver
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
          return new MtlFetiPrimalSolver(this->name);
        }
      };

      /// Constructor of FETI-DP solver class.
      explicit MtlFetiPrimalSolver(string name)
        : LinearSolverInterface(name)
      {}

    };


  }
} // namespace AMDiS

#endif // AMDIS_MTL_FETI_PRIMAL_SOLVER_H
