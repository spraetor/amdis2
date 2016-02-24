/** \file MatrixVector.h */

#pragma once

#include <MatrixVector_fwd.h>
#include <matrix_vector/MemoryBase.hpp>
#include <matrix_vector/Vector.hpp>
#include <matrix_vector/Matrix.hpp>
#include <matrix_vector/ExprConcepts.hpp>

#include <traits/basic.hpp>
#include <Math.hpp>

namespace AMDiS
{
  // ----- traits types --------------------------------------------------------

  namespace detail
  {
    template <class E>
    struct assign_type<E, Requires_t<concepts::Expr<E>>>
    {
      using type = if_then_else<
	/*if */ traits::is_vector<E>::value,
	/* then */ Vector<Value_t<E>>,
	/* else */ if_then_else<
	    /* if */ traits::is_matrix<E>::value,
	    /* then */ Matrix<Value_t<E>>,
	    /* else */ E>>;
    };

  } // end namespace detail


  namespace detail
  {
    // ----- Vector-to-Matrix transform ----------------------------------------

    template <class V>
    struct VectorToMatrix
    {
      STATIC_ASSERT( traits::is_vector<V>::value );

      using type = Matrix<Value_t<V>>;
    };

    template <class T, small_t N>
    struct VectorToMatrix<StaticVector<T,N>>
    {
      using type = StaticMatrix<T, N, N>;
    };

    template <class T>
    struct VectorToMatrix<DimVec<T>>
    {
      using type = DimMat<T>;
    };

    template <class T>
    struct VectorToMatrix<WorldVector<T>>
    {
      using type = WorldMatrix<T>;
    };

    // ----- Matrix-to-Vector transform ----------------------------------------

    template <class M>
    struct MatrixToVector
    {
      STATIC_ASSERT( traits::is_matrix<M>::value );

      using type = Vector<Value_t<M>>;
    };

    template <class T, small_t N>
    struct MatrixToVector<StaticMatrix<T,N,N>>
    {
      using type = StaticVector<T, N>;
    };

    template <class T>
    struct MatrixToVector<DimMat<T>>
    {
      using type = DimVec<T>;
    };

    template <class T>
    struct MatrixToVector<WorldMatrix<T>>
    {
      using type = WorldVector<T>;
    };
  }

  template <class V>
  using VectorToMatrix_t = typename detail::VectorToMatrix<V>::type;

  template <class M>
  using MatrixToVector_t = typename detail::MatrixToVector<M>::type;

} // end namespace AMDiS
