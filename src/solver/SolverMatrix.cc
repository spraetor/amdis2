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


#include "solver/SolverMatrix.h"
#include "solver/Mapper.h"
#include "MatrixStreams.h"

namespace AMDiS
{

  void SolverMatrix<Matrix<DOFMatrix*>>::buildMatrix() const
  {
    BlockMapper mapper(*this);
    matrix.change_dim(mapper.getNumRows(), mapper.getNumCols());
    set_to_zero(matrix);
    MatMap<const SolverMatrix<Matrix<DOFMatrix*>>, BlockMapper> matMap(*this,mapper);
    matrix << matMap;

  }

} // end namespace AMDiS
