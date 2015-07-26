/** \file basic.hpp */

#pragma once

#include <boost/version.hpp>

#include "boost/numeric/ublas/detail/returntype_deduction.hpp"
#if BOOST_VERSION >= 105600
#include <boost/core/enable_if.hpp>
#else
#include <boost/utility/enable_if.hpp>
#endif

#include <utility>
#include <type_traits>

#include "meta_basic.hpp"

namespace AMDiS 
{
  
  template <class T>
  using Value_t = typename T::value_type;
  
  template <class T>
  using Size_t = typename T::size_type;
  
  template <class T>
  using Result_t = typename T::result_type;

  // introduce some shortcuts for boost::mpl
  // ---------------------------------------
  using boost::enable_if;
  using boost::enable_if_c;
  using boost::disable_if;
  using boost::disable_if_c;
  
#if 0
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
#endif
  
  
  namespace traits 
  {
  
    // dummy type
    typedef boost::numeric::ublas::error_cant_deduce_type no_valid_type;
  
    template <class A, class B>
    struct is_multiplicable : not_< 
      std::is_same< typename mtl::Multiplicable<A,B>::result_type, 
		                no_valid_type > > {};
  
    template <class A, class B>
    struct is_addable : not_< 
      std::is_same< typename mtl::Addable<A,B>::result_type, 
			              no_valid_type > > {};
      
    template <class T>
    struct is_trivially_copyable : std::is_trivially_copyable<T> {};

  } // end namespace AMDiS
  
} // end namespace AMDiS
