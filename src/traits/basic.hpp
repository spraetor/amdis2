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
  -> decltype(__VA_ARGS__) { return (__VA_ARGS__); }

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

  template <class T1, class T2>
  using Common_t = typename std::common_type<T1,T2>::type;
  
  namespace detail
  {
    template <class T, class = void>
    struct assign_type { using type = T; };
    
  } // end namespace detail
  
  template <class T>
  using Assign_t = typename detail::assign_type<T>::type;
  
  
  // ---------------------------------------------------------------------------
  
  
  template <class... Ts>
  struct Types {};
  
  template <class... Ts>
  using Types_t = Types<Decay_t<Ts>...>;
  
  template <int... Is>
  struct Ints {};
  
  
  // ---------------------------------------------------------------------------
  
  
  using boost::enable_if;
  using boost::enable_if_c;
  using boost::disable_if;
  using boost::disable_if_c;
  
  template <class C, class T = void>
  using Requires_t = typename std::enable_if<C::value,T>::type;
  
  template <bool C, class T = void>
  using Requires = typename std::enable_if<C,T>::type;
  
  // ---------------------------------------------------------------------------
  
  
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
            
    namespace detail
    {
      template <class A, class B>
      struct IsCompatible : std::is_same<Decay_t<A>, Decay_t<B> > {};
      
      template <class A, class B>
      struct IsCompatible<Types<A>, Types<B>> 
        : IsCompatible<A,B> {};
        
      template <>
      struct IsCompatible<Types<>, Types<>> : true_ {};
      
      template <class A0, class... As, class B0, class... Bs>
      struct IsCompatible<Types<A0,As...>, Types<B0,Bs...>>
        : and_<IsCompatible<A0,B0>, IsCompatible<Types<As...>, Types<Bs...>> > {};
    }
    
    template <class A, class B>
    using IsCompatible = detail::IsCompatible<A,B>;
    
  } // end namespace traits
  
  
} // end namespace AMDiS
