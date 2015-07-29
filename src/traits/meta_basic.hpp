/** \file meta_basic.hpp */

#pragma once

#include <type_traits>

#define STATIC_ASSERT(...) \
  static_assert(__VA_ARGS__, #__VA_ARGS__)

namespace AMDiS 
{
  // introduce some shortcuts for integral constants
  // -----------------------------------------------
  template <int I> 
  using int_ = std::integral_constant<int, I>;
  
  static constexpr int_<0> _0 {};
  static constexpr int_<1> _1 {};
  static constexpr int_<2> _2 {};
  static constexpr int_<3> _3 {};
  
  template <bool B> 
  using bool_ = std::integral_constant<bool, B>;
  
  using true_ = bool_<true>;
  using false_ = bool_<false>;
    
  // --- some boolean operations -----------------------------------------------
  namespace detail
  {
    template <class T0, class... Ts>
    struct or_ : bool_<T0::value || or_<Ts...>::value> {};
    
    template <class T0>
    struct or_<T0> : bool_<T0::value> {};
      
    template <class T0, class... Ts>
    struct and_ : bool_<T0::value && and_<Ts...>::value> {};
    
    template <class T0>
    struct and_<T0> : bool_<T0::value> {};
    
  } // end namespace detail
  
  template <class... Ts>
  using and_ = detail::and_<Ts...>;
  
  template <bool... Bs>
  using and_c = and_<bool_<Bs>...>;
  
  //template <class... Ts>
  //static constexpr bool and_v = and_<Ts...>::value;
  
  template <class... Ts>
  using or_ = detail::or_<Ts...>;
  
  template <bool... Bs>
  using or_c = or_<bool_<Bs>...>;
  
  //template <class... Ts>
  //static constexpr bool or_v = or_<Ts...>::value;
  
  template <class A>
  using not_ = bool_<!(A::value)>;
  
  //template <class T>
  //static constexpr bool not_v = not_<T>::value;
  
  
  template <bool C, class T1, class T2>
  using if_then_else = typename std::conditional<C, T1, T2>::type;
  
  template <class T, T A, T B>
  using is_equal = if_then_else< A == B, true_, false_ >;
  
  template <class T, class... Ts>
  using is_one_of = or_< std::is_same<T, Ts>... >;
  
  
  /// for_each for std::tuple
  template <int I = 0, class FuncT, class... Tp>
  inline typename std::enable_if<(I == sizeof...(Tp)), void>::type
  for_each(std::tuple<Tp...> &, FuncT) { }

  template <int I = 0, class FuncT, class... Tp>
  inline typename std::enable_if<(I < sizeof...(Tp)), void>::type
  for_each(std::tuple<Tp...>& t, FuncT f)
  {
    f(std::get<I>(t));
    for_each< (I + 1), FuncT, Tp...>(t, f);
  }
  
  
  namespace detail
  {
    template <class T, class = decltype(&T::operator[])>
    static true_ isVectorImpl(int);
    
    template <typename T>
    static false_ isVectorImpl(...);
    
  } // end namespace detail
  
  template <class T>
  struct IsVector
    : decltype(detail::isVectorImpl<T>(int{}))
  {};
} // end namespace AMDiS
