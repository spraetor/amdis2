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



/** \file expr_traits.hpp */

#ifndef AMDIS_EXPR_TRAITS_HPP
#define AMDIS_EXPR_TRAITS_HPP

#include "Traits.h"
#include "LazyOperatorTerm.h"
#include "value_expr.hpp"

namespace AMDiS 
{
  namespace traits
  {    
      
    // type-trait for short-cuts for constants
    // ___________________________________________________________________________  
	  
    template<typename T>
    struct is_constant : is_numeric<T>::type {};
    
    template<typename T>
    struct is_constant<WorldVector<T> > : boost::mpl::bool_<true> {};
    
    template<typename T>
    struct is_constant<WorldMatrix<T> > : boost::mpl::bool_<true> {};
    
    
    // type-traits for terms
    // ___________________________________________________________________________
    template<typename T>
    struct is_expr : boost::is_base_of<LazyOperatorTermBase, T>::type {};
    
    
    // type-traits for arguments, filter terms and constants
    // ___________________________________________________________________________
    template<typename T>
    struct is_valid_arg : boost::mpl::or_
    <
      typename is_expr<T>::type,
      typename is_constant<T>::type
    >::type {};
    
    template<typename T1, typename T2>
    struct is_valid_arg2 : boost::mpl::and_
    <
      typename is_valid_arg<T1>::type,
      typename is_valid_arg<T2>::type,
      typename boost::mpl::or_
      <
	typename is_expr<T1>::type,
	typename is_expr<T2>::type
      >::type
    >::type {};
    
    template<typename T1, typename T2, typename T3>
    struct is_valid_arg3 : boost::mpl::and_
    <
      typename is_valid_arg<T1>::type,
      typename is_valid_arg<T2>::type,
      typename is_valid_arg<T3>::type,
      typename boost::mpl::or_
      <
	typename is_expr<T1>::type,
	typename is_expr<T2>::type,
	typename is_expr<T3>::type
      >::type
    >::type {};
    
    
      
    // expressions
    template < typename T >
    struct category<T, typename boost::enable_if< typename is_expr<T>::type >::type >
    {
      typedef tag::expression         tag;
      typedef typename T::value_type  value_type;
      // typedef size_t               size_type;
    };
	
      
    // type-conversion
    // ___________________________________________________________________________
    
    template<typename T, bool IsValid, bool IsConstant>
    struct to_expr_aux {
      typedef void type;
    };
      
    template<typename T>
    struct to_expr_aux<T, /*IsValid: */true, /*IsConstant: */false> 
    {
      typedef T type;
      static const type& get(const T& t)
      {
	return t;
      }
    };
    
    template<typename T> // T is numeric constant
    struct to_expr_aux<T, /*IsValid: */true, /*IsConstant: */true> 
    {
      typedef ::AMDiS::expressions::RValue<T> type;
      static type get(const T& t)
      {
	return type(t);
      }
    };
      
    template<typename T>
    struct to_expr 
    {
      typedef to_expr_aux<T, is_valid_arg<T>::value, is_constant<T>::value> to;
      typedef typename to::type type;
    };
    
  } // end namespace traits

} // end namespace AMDiS

#endif // AMDIS_EXPR_TRAITS_HPP
