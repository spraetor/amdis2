/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file Traits.h */

#ifndef AMDIS_TRAITS_H
#define AMDIS_TRAITS_H

#include <iostream>
#include "FixVec.h"
#include "MatrixVector.h"
#include "AMDiS_fwd.h"

#include "traits/tag.hpp"
#include "traits/types.hpp"
#include "traits/category.hpp"
#include "traits/num_rows.hpp"
#include "traits/num_cols.hpp"
#include "traits/size.hpp"

#include "boost/numeric/ublas/detail/returntype_deduction.hpp"
#include "boost/numeric/mtl/concept/std_concept.hpp" // for return_type deduction

// TODO: move to traits directory

namespace AMDiS 
{
  namespace traits
  {
    
    // type-generators
    // _________________________________________________________________________
    
    // larger types
    template<typename T1, typename T2>
    struct larger_type
    {
      typedef typename boost::mpl::if_c< (sizeof(T1) > sizeof(T2)), T1, T2 >::type type;
    };
    
    // TODO: improve this implementation!!
#if HAS_DECLTYPE
    
    /// determines the type of the product T1*T2
    template<typename T1, typename T2>
    struct mult_type
    {      
      typedef boost::numeric::ublas::type_deduction_detail::base_result_of<T1, T2> base_type;
      static typename base_type::x_type x;
      static typename base_type::y_type y;
      
      typedef decltype( x * y ) type;
    };
    
    /// determines the type of the sum T1+T2
    template<typename T1, typename T2>
    struct add_type
    {
      typedef boost::numeric::ublas::type_deduction_detail::base_result_of<T1, T2> base_type;
      static typename base_type::x_type x;
      static typename base_type::y_type y;
      
      typedef decltype( x + y ) type;
    };
    
#else
    
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
      typedef typename mtl::Multiplicable<T1, T2>::result_type type;
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
      typedef typename mtl::Addable<T1, T2>::result_type type;
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
    
#endif // endif(HAS_DECLTYPE)
      
      
  } // end namespace traits

} // end namespace AMDiS

#endif
