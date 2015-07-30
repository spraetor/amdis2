/** \file Traits.h */

#pragma once

// std c++ headers
#include <iostream>

// AMDiS headers
#include "AMDiS_fwd.h"

#include "traits/basic.hpp"
#include "traits/meta_basic.hpp"
#include "traits/tag.hpp"
#include "traits/category.hpp"


namespace AMDiS 
{
  namespace traits
  {        
    namespace detail
    {
      template <bool Valid, class Type>
      struct Enabler
      {
        using type = typename Type::type;
      };
      
      template <class Type>
      struct Enabler<false, Type>
      {
        using type = no_valid_type;
      };
      
    } // end namespace detail
    
    // multiplication types
    // _________________________________________________________________________
    template<typename T1, typename T2, typename Category1, typename Category2>
    struct mult_type_dispatch
    {
      typedef no_valid_type type;
    };
    
    /// determines the type of the product T1*T2
    template<typename T1, typename T2>
    struct mult_type : mult_type_dispatch
      <
      	T1, T2,
      	typename category<T1>::tag, 
      	typename category<T2>::tag
      > {};
    
    /// Scalar*Scalar => Scalar
    template<typename T1, typename T2>
    struct mult_type_dispatch<T1, T2, tag::scalar, tag::scalar>
    {
      template <class T1_, class T2_>
      struct Type { using type = decltype(std::declval<T1_>() * std::declval<T2_>()); };
      using type = typename detail::Enabler< is_multiplicable<T1,T2>::value, Type<T1,T2> >::type;
    };
    
    /// Vec*Vec => Scalar (dot-product)
    template<typename T1, typename T2>
    struct mult_type_dispatch<T1, T2, tag::vector, tag::vector>
    {
      typedef typename mult_type
      <
      	typename category<T1>::value_type,
      	typename category<T2>::value_type
      >::type type;
    };
    
    /// Mat*Mat => Mat
    template<typename T>
    struct mult_type_dispatch<T, T, tag::matrix, tag::matrix>
    {
      typedef T type;
    };
    
    
    /// Vec*Scalar => Vector
    template<template<class> class Container, typename T1, typename T2>
    struct mult_type_dispatch<Container<T1>, T2, tag::vector, tag::scalar>
    {
      typedef typename mult_type<T1, T2>::type value_type;
      typedef Container<value_type> type;
    };
    
    /// Scalar*Vector => Vector
    template<template<class> class Container, typename T1, typename T2>
    struct mult_type_dispatch<T1, Container<T2>, tag::scalar, tag::vector>
    {
      typedef typename mult_type<T1, T2>::type value_type;
      typedef Container<value_type> type;
    };
    
    
    /// Matrix*Scalar => Matrix
    template<template<class> class Container, typename T1, typename T2>
    struct mult_type_dispatch<Container<T1>, T2, tag::matrix, tag::scalar>
    {
      typedef typename mult_type<T1, T2>::type value_type;
      typedef Container<value_type> type;
    };
    
    /// Scalar*Matrix => Matrix
    template<template<class> class Container, typename T1, typename T2>
    struct mult_type_dispatch<T1, Container<T2>, tag::scalar, tag::matrix>
    {
      typedef typename mult_type<T1, T2>::type value_type;
      typedef Container<value_type> type;
    };
    
    
    /// Matrix*Vector => Vector
    template<template<class> class MatrixType, typename T1, template<class> class VectorType, typename T2>
    struct mult_type_dispatch<MatrixType<T1>, VectorType<T2>, tag::matrix, tag::vector>
    {
      typedef typename mult_type<T1, T2>::type value_type;
      typedef VectorType<value_type> type;
    };
    
    
    // addition types
    // _________________________________________________________________________
    template<typename T1, typename T2, typename Category1, typename Category2>
    struct add_type_dispatch
    {
      typedef no_valid_type type;
    };
    
    /// determines the type of the sum T1+T2
    template<typename T1, typename T2>
    struct add_type : add_type_dispatch
      <
      	T1, T2,
      	typename category<T1>::tag, 
      	typename category<T2>::tag
      > {};
    
    /// Scalar+Scalar => Scalar
    template<typename T1, typename T2>
    struct add_type_dispatch<T1, T2, tag::scalar, tag::scalar>
    {
      template <class T1_, class T2_>
      struct Type { using type = decltype(std::declval<T1_>() + std::declval<T2_>()); };
      using type = typename detail::Enabler< is_addable<T1,T2>::value, Type<T1,T2> >::type;
    };
    
    /// Vec+Vec => Vec
    template<template<class> class Container, typename T1, typename T2>
    struct add_type_dispatch<Container<T1>, Container<T2>, tag::vector, tag::vector>
    {
      typedef typename add_type<T1, T2>::type value_type;
      typedef Container<value_type> type;
    };
    
    /// Mat+Mat => Mat
    template<template<class> class Container, typename T1, typename T2>
    struct add_type_dispatch<Container<T1>, Container<T2>, tag::matrix, tag::matrix>
    {
      typedef typename add_type<T1, T2>::type value_type;
      typedef Container<value_type> type;
    };
    
    
    /// Vec+Scalar => Vector
    template<template<class> class Container, typename T1, typename T2>
    struct add_type_dispatch<Container<T1>, T2, tag::vector, tag::scalar>
    {
      typedef typename add_type<T1, T2>::type value_type;
      typedef Container<value_type> type;
    };    
    
    /// Mat+Scalar => Mat
    template<template<class> class Container, typename T1, typename T2>
    struct add_type_dispatch<Container<T1>, T2, tag::matrix, tag::scalar>
    {
      typedef typename add_type<T1, T2>::type value_type;
      typedef Container<value_type> type;
    };
    
      
  } // end namespace traits

} // end namespace AMDiS
