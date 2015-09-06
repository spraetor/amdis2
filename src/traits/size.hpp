/** \file size.hpp */

#pragma once

#include <AMDiS_fwd.h>
#include <traits/traits_fwd.hpp>
#include <traits/scalar_types.hpp>

namespace AMDiS
{
  template <class T, class = Requires_t<concepts::Arithmetic<T>>>
  inline size_t size(T const& v)
  {
    return 1;
  }
  
  template <class FixVecType>
  inline size_t size(VectorOfFixVecs<FixVecType> const& v)
  {
    return v.getSize() * v.getSizeOfFixVec();
  }
  
  template <class FixVecType>
  inline size_t size(MatrixOfFixVecs<FixVecType> const& v)
  {
    return v.getNumberOfRows() == 0 
      ? 0 
      : v.getNumberOfRows() * v.getNumberOfColumns() * v[0].getSizeOfFixVec();
  }

} // end namespace AMDiS
