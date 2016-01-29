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



/** \file BddcMlSolver.h */

#ifndef AMDIS_BDDC_ML_SOLVER_H
#define AMDIS_BDDC_ML_SOLVER_H

#include "AMDiS_fwd.h"
#include "parallel/ParallelSolver.h"

namespace AMDiS
{
  namespace Parallel
  {

#ifdef HAVE_BDDC_ML

    class BddcMlSolver : public ParallelSolver
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
          return new BddcMlSolver(this->name);
        }
      };

      BddcMlSolver(std::string name)
        : ParallelSolver(name, true) {}

    protected:
      /// Implementation of \ref LinearSolverInterface::solveSystemImpl()
      virtual int solveSystemImpl(SolverMatrix<Matrix<DOFMatrix*>> const& A,
				  SystemVector& x,
				  SystemVector& b,
				  bool createMatrixData,
				  bool storeMatrixData) override;

      void addDofMatrix(const DOFMatrix* mat,
                        std::vector<int>& i_sparse,
                        std::vector<int>& j_sparse,
                        std::vector<double>& a_sparse,
                        int nComponents,
                        int ithRowComponent,
                        int ithColComponent);
    };

#endif

  } // end namespace Parallel
} // end namespace AMDiS

#endif

