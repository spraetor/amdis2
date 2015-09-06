/** \file category.hpp */

#pragma once

#include <type_traits>

#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/coordinate2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp>

// AMDiS headers
#include <AMDiS_fwd.h>
#include <AMDiS_base.h>
#include <traits/traits_fwd.hpp>
#include "tag.hpp"

namespace AMDiS
{
  namespace traits
  {
    // scalars
    template <class T>
    struct category<T, Requires_t<concepts::Arithmetic<T>> >
    {
      typedef tag::scalar tag;
      typedef T           value_type;
      //       typedef size_t      size_type;
    };

    
    // vectors
//     template <class T>
//     struct category<T, Requires_t<and_<traits::HasValueType<T>, traits::HasSizeType<T>, concepts::Vector<T>>> >
//     {
//       using tag = tag::vector;
//       using value_type = Value_t<T>;
//       using size_type  = Size_t<T>;
//     };
// 
//     
//     // matrices
//     template <class T>
//     struct category<T, Requires_t<and_<traits::HasValueType<T>, traits::HasSizeType<T>, concepts::Matrix<T>>> >
//     {
//       using tag = tag::matrix;
//       using value_type = Value_t<T>;
//       using size_type  = Size_t<T>;
//     };

    template <>
    struct category<AMDiS::SystemVector>
    {
      typedef tag::vector      tag;
      typedef double           value_type;
      typedef int              size_type;
    };

    template <class T>
    struct category<WorldVector<WorldVector<T>>>
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
      typedef int          size_type;
    };

    template <class Matrix>
    struct category<AMDiS::SolverMatrix<Matrix>>
    {
      typedef tag::matrix  tag;
    };


    template <class T>
    struct category<const T> : category<T> {};

    // -------------------------------------------------------------------------

    template <class T, class S, class = void>
    struct IsConvertible : std::is_convertible<T, S> {};

    template <class T, class S>
    struct IsConvertible<T, S, Requires_t<and_<is_vector<T>, is_vector<S>>>>
      : std::is_convertible<Value_t<T>, Value_t<S>> {};

    template <class T, class S>
    struct IsConvertible<T, S, Requires_t<and_<is_matrix<T>, is_matrix<S>>>>
      : std::is_convertible<Value_t<T>, Value_t<S>> {};

  } // end namespace traits

} // end namespace AMDiS
