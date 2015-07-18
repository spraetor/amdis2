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



/** \file basic.hpp */

#ifndef AMDIS_BASIC_TRAITS_HPP
#define AMDIS_BASIC_TRAITS_HPP

#include <boost/version.hpp> 
#include <boost/mpl/logical.hpp>

#include "boost/numeric/ublas/detail/returntype_deduction.hpp"
#if BOOST_VERSION >= 105600
#include <boost/core/enable_if.hpp>
#else
#include <boost/utility/enable_if.hpp>
#endif

#include <utility>
#include <type_traits>

namespace AMDiS 
{

  // introduce some shortcuts for boost::mpl
  // ---------------------------------------
  using boost::mpl::bool_;
  using boost::mpl::true_;
  using boost::mpl::false_;
  using boost::mpl::and_;
  using boost::mpl::or_;
  
  using boost::enable_if;
  using boost::enable_if_c;
  using boost::disable_if;
  using boost::disable_if_c;
  
  template <class... Iters>
  struct MultiIterator
  {
    template <class... Iters_>
    MultiIterator(Iters_&&... iters_) 
      : iters(std::forward<Iters_>(iters_)...) { }
    
    auto begin() {
      using Indices = std::make_index_sequence<sizeof...(Iters)>;
      return begin_impl(Indices());
    }
    
    auto end() {
      using Indices = std::make_index_sequence<sizeof...(Iters)>;
      return end_impl(Indices());
    }
    
  private:    
    template <size_t ... I>
    auto begin_impl(std::index_sequence<I...>) {
      return std::make_tuple(std::get<I>(iters).begin()...);
    }
    template <size_t ... I>
    auto end_impl(std::index_sequence<I...>) {
      return std::make_tuple(std::get<I>(iters).end()...);
    }
    
    using TupleType = std::tuple<Iters...>;
    TupleType iters;
  };
  
  
  template <class... Iters>
  auto make_iter(Iters&&... iters)
  {
    return MultiIterator<Iters...>(std::forward<Iters>(iters)...);
  }
  
  
  
  namespace traits 
  {
  
    // dummy type
    typedef boost::numeric::ublas::error_cant_deduce_type no_valid_type;
  
    template <class A, class B>
    struct is_multiplicable : boost::mpl::not_< 
	boost::is_same< typename mtl::Multiplicable<A,B>::result_type, 
		        no_valid_type > > {};
  
    template <class A, class B>
    struct is_addable : boost::mpl::not_< 
	boost::is_same< typename mtl::Addable<A,B>::result_type, 
			no_valid_type > > {};
      
#ifdef HAS_CPP11
    template <typename T>
    struct is_trivially_copyable : std::is_trivially_copyable<T> {};
#else
    template <typename T>
    struct is_trivially_copyable : boost::is_pod<T> {};
#endif
    
    template <class T, T A, T B>
    struct equal : boost::mpl::if_c< A == B, true_, false_ > {};
  }
  
//   template <class T> 
//   T zero(T result) { result = 0; return result; }
}


#endif // AMDIS_BASIC_TRAITS_HPP
