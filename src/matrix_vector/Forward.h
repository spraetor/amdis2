/** \file Forward.h */

#pragma once

#include <Config.h>		// small_t

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

namespace AMDiS
{
  // some forwards declaration

  // memory-policies
  template <class T, small_t N, small_t M>   struct MemoryBaseStatic;
  template <class T, bool aligned>          struct MemoryBaseDynamic;
  template <class T, small_t N, small_t M>   struct MemoryBaseHybrid;

  // size-policies
  struct DefaultSizePolicy;
  template <size_t S> struct StaticSizePolicy;

  // matrix-vector types
  template <class Model, class MemoryPolicy>   struct MatrixVectorBase;
  template <class MemoryPolicy, class SizePolicy>    struct VectorBase;
  template <class MemoryPolicy, class SizePolicy>    struct MatrixBase;


  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    template <class M, class S>
    struct category<VectorBase<M,S>>
    {
      typedef tag::vector        tag;
      typedef Value_t<M>  value_type;
      typedef Size_t<M>    size_type;
    };

    template <class M, class S>
    struct category<MatrixBase<M,S>>
    {
      typedef tag::matrix        tag;
      typedef Value_t<M>  value_type;
      typedef Size_t<M>    size_type;
    };
    /// \endcond

  } // end namespace traits
} // end namespace AMDiS
