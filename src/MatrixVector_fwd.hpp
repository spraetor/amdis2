#pragma once

#include "GeoIndex.hpp"
#include "matrix_vector/Forward.hpp"

namespace AMDiS
{
  // size-policy forward declaration
  template <GeoIndex G> struct FixedSize;
  template <GeoIndex G> struct MaxSize;  // maybe a definition for template parameters must be available here


  // ----- Vector types --------------------------------------------------------

  /// define a FixVec as a specialized static-vector
#if AMDIS_FIXED_SIZE
  template <class T, GeoIndex G> 
  using FixVec
    = VectorBase<MemoryBaseHybrid<T, MaxSize<G>::value, 1>, FixedSize<G>>;
#else
  template <class T, GeoIndex G> 
  using FixVec
    = VectorBase<MemoryBaseDynamic<T, false>, FixedSize<G>>;
#endif

  /// define a WorldVector as a specialized FixVec
  template <class T> 
  using WorldVector = FixVec<T, WORLD>;

  /// define a DimVec as a specialized FixVec
  template <class T> 
  using DimVec = FixVec<T, PARTS>;

  template <class T, small_t N>
  using StaticVector = VectorBase<MemoryBaseStatic<T, N, 1>, StaticSizePolicy<N>>;

  /// define a Vector as a specialized dynamic-vector
  template <class T> 
  using Vector
    = VectorBase<MemoryBaseDynamic<T, false>, DefaultSizePolicy>;

  // ----- Matrix types --------------------------------------------------------

  /// define a FixMat as a specialized static-matrix
#if AMDIS_FIXED_SIZE
  template <class T, GeoIndex G> 
  using FixMat
    = MatrixBase<MemoryBaseHybrid<T, MaxSize<G>::value, MaxSize<G>::value>, FixedSize<G>>;
#else
  template <class T, GeoIndex G> 
  using FixMat
    = MatrixBase<MemoryBaseDynamic<T, false>, FixedSize<G>>;
#endif

  /// define a WorldMatrix as a specialized FixMat
  template <class T> 
  using WorldMatrix = FixMat<T, WORLD>;

  /// define a DimMat as a specialized FixMat
  template <class T> 
  using DimMat = FixMat<T, PARTS>;

  template <class T, small_t N, small_t M>
  using StaticMatrix = MatrixBase<MemoryBaseStatic<T, N, M>, DefaultSizePolicy>;

  /// define a Matrix as a specialized dynamic-matrix
  template <class T> 
  using Matrix
    = MatrixBase<MemoryBaseDynamic<T, false>, DefaultSizePolicy>;

} // end namespace AMDiS
