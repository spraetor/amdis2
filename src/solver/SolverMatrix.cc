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
    MatMap<const SolverMatrix<Matrix<DOFMatrix*>>, BlockMapper> matMap{*this, mapper};
    matrix << matMap;

  }

} // end namespace AMDiS
