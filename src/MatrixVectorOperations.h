/** \file MatrixVectorOperations.h */

#pragma once

#include <Math.h>
#include <MatrixVector.h>
#include <matrix_vector/MatrixVectorOperations.hpp>

namespace AMDiS
{
  template <class T> inline void nullify(WorldVector<T>& a)
  {
    T null;
    nullify(null);
    a = null;
  }

  template <class T> inline void nullify(WorldMatrix<T>& a)
  {
    T null;
    nullify(null);
    a = null;
  }

} // end namespace AMDiS
