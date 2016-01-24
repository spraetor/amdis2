/** \file size.hpp */

#pragma once

#include <AMDiS_fwd.h>
#include <traits/traits_fwd.hpp>
#include <traits/scalar_types.hpp>

#include <tuple>
#include <array>

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
  
  
  // Size-function for compile-time containers
    
  template <class Tuple>
  struct Size : int_<std::tuple_size<Tuple>::value> {};
  
  template <template<class...> class Tuple, class... Ts>
  struct Size<Tuple<Ts...>> : int_<sizeof...(Ts)> {};
  
  template <class T, size_t N>
  struct Size<std::array<T, N>> : int_<N> {};

} // end namespace AMDiS
