#pragma once

// AMDiS includes
#include "Config.hpp"		// small_t

#include "traits/basic.hpp"
#include "traits/traits_fwd.hpp"

namespace AMDiS
{
  // some forwards declaration

  // memory-policies
  template <class T, small_t N, small_t M>  struct MemoryBaseStatic;
  template <class T, bool aligned>          struct MemoryBaseDynamic;
  template <class T, small_t N, small_t M>  struct MemoryBaseHybrid;

  // size-policies
  struct DefaultSizePolicy;
  template <size_t S> struct StaticSizePolicy;

  // matrix-vector types
  template <class Model, class MemoryPolicy>       struct MatrixVectorBase;
  template <class MemoryPolicy, class SizePolicy>  struct VectorBase;
//   template <class MemoryPolicy, class SizePolicy>  struct MatrixBase;


  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    template <class M, class S>
    struct category<VectorBase<M,S>>
    {
      using tag = tag::vector;
      using value_type = Value_t<M>;
      using size_type  = Size_t<M>;
    };

//     template <class M, class S>
//     struct category<MatrixBase<M,S>>
//     {
//       using tag = tag::matrix;
//       using value_type = Value_t<M>;
//       using size_type  = Size_t<M>;
//     };
    /// \endcond

  } // end namespace traits
} // end namespace AMDiS
