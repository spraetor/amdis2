/** \file basic.hpp */

#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// boost headers
#include <boost/utility/enable_if.hpp>

// AMDiS headers
#include "meta_basic.hpp"


// a helper mycro to reduce typing
#define RETURNS(...) \
  noexcept(noexcept(__VA_ARGS__)) \
    -> decltype(__VA_ARGS__) { return (__VA_ARGS__); }
    
#define RETURNS_CONST(...) \
  noexcept(noexcept(__VA_ARGS__)) \
    -> decltype(__VA_ARGS__) const { return (__VA_ARGS__); }

namespace AMDiS 
{
  template <class T>
  using Value_t = typename T::value_type;
  
  template <class T>
  using Size_t = typename T::size_type;
  
  template <class T>
  using Result_t = typename T::result_type;
  
  template <class T>
  using Decay_t = typename std::decay<T>::type;

  // introduce some shortcuts for boost::mpl
  // ---------------------------------------
  using boost::enable_if;
  using boost::enable_if_c;
  using boost::disable_if;
  using boost::disable_if_c;
  
  // some traits to test for binary operators on types
  namespace traits 
  {
    template <class T>
    struct id_
    {
      using type = T;
    };
    
    namespace detail
    {
      template <class...>
      struct voider { using type = void; };
      
      template <class...Ts>
      using void_t = typename voider<Ts...>::type;
  
      template <class T, class U, class BinOp, class Enable = void>
      struct IsBinopAbleImpl : false_ {};
      
      template <class T, class U, class BinOp>
      struct IsBinopAbleImpl<T, U, BinOp, void_t<
        decltype(std::declval<BinOp>()(std::declval<T>(), 
                                       std::declval<U>() )) > > : true_ {};

    } // end namespace detail
    
    template <class T, class U, class BinOp>
    using IsBinopAble = typename detail::IsBinopAbleImpl<T, U, BinOp>::type;
  
    namespace detail
    {
      struct PlusOp {
        template <class T, class U>
        auto operator()(T t, U u) RETURNS(t + u)
      };
      
      struct MultipliesOp {
        template <class T, class U>
        auto operator()(T t, U u) RETURNS(t * u)
      };

    } // end namespace detail
    
    template <class T, class U = T>
    using is_addable = IsBinopAble<T, U, detail::PlusOp>;
    
    template <class T, class U = T>
    using is_multiplicable = IsBinopAble<T, U, detail::MultipliesOp>;
      
    template <class T>
    using is_trivially_copyable = std::is_pod<T>;

    // larger types
    template <class... Ts>
    struct larger_type;
    
    template <class T1, class T2, class... Ts>
    struct larger_type<T1, T2, Ts...>
    {
      using type = typename if_then_else< (sizeof(T1) > sizeof(T2)), 
					  larger_type<T1,Ts...>, 
					  larger_type<T2,Ts...> >::type;
    };
    
    template <class T1, class T2>
    struct larger_type<T1, T2>
    {
      using type = if_then_else<(sizeof(T1) > sizeof(T2)), T1, T2>;
    };
    
    template <class T> struct larger_type<T> { using type = T; };
    
    // maximal size type
    template <class... Es>
    using max_size_type = typename larger_type<Size_t<Es>...>::type;
    
  } // end namespace traits
  
} // end namespace AMDiS
