#pragma once

#include "Math.hpp"
#include "MatrixVector.hpp"
#include "matrix_vector/MatrixVectorOperations.hpp"

namespace AMDiS
{
  template <class T> inline void nullify(WorldVector<T>& a)
  {
    T null;
    nullify(null);
    a = null;
  }

} // end namespace AMDiS
