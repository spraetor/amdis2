/** \file meta_basic.hpp */

#pragma once

// std c++ headers
#include <type_traits>

namespace AMDiS
{
  // introduce some shortcuts for integral constants
  // ---------------------------------------------------------------------------

  template <int I>
  using int_ = std::integral_constant<int, I>;

  static constexpr int_<0> _0 {};
  static constexpr int_<1> _1 {};
  static constexpr int_<2> _2 {};
  static constexpr int_<3> _3 {};

  template <bool B>
  using bool_  = std::integral_constant<bool, B>;
  using true_  = bool_<true>;
  using false_ = bool_<false>;

  static constexpr bool_<true>  _true  {};
  static constexpr bool_<false> _false {};

  template <int I, int J>
  struct range_
  {
    using type = range_;

    static constexpr int  begin()  { return I; }
    static constexpr int  end()    { return J; }
    static constexpr bool empty()  { return I >= J; }

    constexpr operator bool() const { return I < J; }
  };

  template <int I>
  using empty_range_ = range_<I,I>;


  // some boolean operations
  // ---------------------------------------------------------------------------

  namespace detail
  {
    template <class... Ts> struct or_;
    
    template <class T0, class... Ts>
    struct or_<T0, Ts...> : bool_<T0::value || or_<Ts...>::value> {};

    template <>
    struct or_<> : false_ {};

    template <class... Ts> struct and_;
    
    template <class T0, class... Ts>
    struct and_<T0, Ts...> : bool_<T0::value && and_<Ts...>::value> {};

    template <>
    struct and_<> : true_ {};

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
  using is_equal = if_then_else<A == B, true_, false_>;

  template <class T, class... Ts>
  using is_one_of = or_<std::is_same<T, Ts>...>;

} // end namespace AMDiS
