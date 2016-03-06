#pragma once

#include "GeoIndex.hpp"
#include "matrix_vector/Forward.hpp"

namespace AMDiS
{
  // size-policy forward declaration
  template <GeoIndex G> struct FixedSize;
  template <GeoIndex G> struct MaxSize;  // maybe a definition for template parameters must be available here


  // ----- Vector types --------------------------------------------------------

  /// define a WorldVector as a specialized FixVec
  template <class T>
  using WorldVector = VectorBase<MemoryBaseDynamic<T, false>, FixedSize<WORLD>>;

} // end namespace AMDiS
