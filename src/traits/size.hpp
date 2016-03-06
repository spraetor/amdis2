/** \file size.hpp */

#pragma once

// std c++ headers
#include <tuple>
#include <array>

// AMDiS headers
#include "AMDiS_fwd.hpp"
#include "traits/traits_fwd.hpp"
#include "traits/scalar_types.hpp"

namespace AMDiS
{
  template <class T>
    Requires_t<concepts::Arithmetic<T>, size_t>
  inline size(T const& /*v*/)
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
