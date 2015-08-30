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



/** \file MatrixNnzStructure.h */

#ifndef AMDIS_MATRIX_NNZ_STRUCTURE_H
#define AMDIS_MATRIX_NNZ_STRUCTURE_H

#include "AMDiS_fwd.h"
#include "parallel/ParallelDofMapping.h"
#include "parallel/PeriodicMap.h"

namespace AMDiS
{
  namespace Parallel
  {

    class MatrixNnzStructure
    {
    public:
      MatrixNnzStructure()
        : dnnz(NULL),
          onnz(NULL)
      {}

      ~MatrixNnzStructure();

      void clear();

      void create(Matrix<DOFMatrix*>& mat,
                  ParallelDofMapping& rowDofMap,
                  ParallelDofMapping& colDofMap,
                  PeriodicMap* perMap,
                  const ElementObjectDatabase& elObjDb,
                  bool localMatrix = false);

      void create(Matrix<DOFMatrix*>& mat,
                  ParallelDofMapping& dofMap,
                  PeriodicMap* perMap,
                  const ElementObjectDatabase& elObjDb,
                  bool localMatrix = false)
      {
        create(mat, dofMap, dofMap, perMap, elObjDb, localMatrix);
      }

    protected:
      void create(int nRows0, int nRows1 = -1);

    public:
      // Array containing the number of nonzeros in the various rows of the DIAGONAL
      // portion of the local submatrix (possibly different for each row)
      int* dnnz;

      // Array containing the number of nonzeros in the various rows of the OFF-DIAGONAL
      // portion of the local submatrix (possibly different for each row)
      int* onnz;
    };

  }
}

#endif
