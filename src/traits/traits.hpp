/** \file traits.hpp */

#pragma once

#include <type_traits>
#include <traits/basic.hpp>

namespace AMDiS
{  
  // some traits to test for binary operators on types
  namespace traits
  {
    namespace detail
    {
      template <class...>
      struct Voider
      {
        using type = void;
      };

      template <class...Ts>
      using Voider_t = typename Voider<Ts...>::type;

      template <class T, class U, class BinOp, class = void>
      struct IsBinopAbleImpl : false_ {};

      template <class T, class U, class BinOp>
      struct IsBinopAbleImpl<T, U, BinOp, Voider_t<
          decltype( std::declval<BinOp>()(std::declval<T>(), std::declval<U>()) )>> 
        : true_ {};

    } // end namespace detail

    template <class T, class U, class BinOp>
    using IsBinopAble = typename detail::IsBinopAbleImpl<T, U, BinOp>::type;

    namespace detail
    {
      struct PlusOp
      {
        template <class T, class U>
        auto operator()(T t, U u) const RETURNS(t + u)
      };

      struct MultipliesOp
      {
        template <class T, class U>
        auto operator()(T t, U u) const RETURNS(t * u)
      };

    } // end namespace detail

    template <class T, class U = T>
    using IsAddable = IsBinopAble<T, U, detail::PlusOp>;

    template <class T, class U = T>
    using IsMultiplicable = IsBinopAble<T, U, detail::MultipliesOp>;

    
    template <class T>
    using IsTriviallyCopyable = std::is_pod<T>;

    // -------------------------------------------------------------------------
    
    // larger types
    template <class... Ts>
    struct LargerType;

    template <class T1, class T2, class... Ts>
    struct LargerType<T1, T2, Ts...>
    {
      using type = typename if_then_else< (sizeof(T1) > sizeof(T2)),
            LargerType<T1,Ts...>,
            LargerType<T2,Ts...> >::type;
    };

    template <class T1, class T2>
    struct LargerType<T1, T2>
    {
      using type = if_then_else<(sizeof(T1) > sizeof(T2)), T1, T2>;
    };

    template <class T> struct LargerType<T>
    {
      using type = T;
    };

    // maximal size type
    template <class... Es>
    using MaxSizeType = typename LargerType<Size_t<Es>...>::type;

    // -------------------------------------------------------------------------
    
    namespace detail
    {
      template <class A, class B>
      struct is_compatible 
        : std::is_same<Decay_t<A>, Decay_t<B>> {};

      template <class A, class B>
      struct is_compatible<Types<A>, Types<B>>
        : is_compatible<A,B> {};

      template <>
      struct is_compatible<Types<>, Types<>> : true_ {};

      template <class A0, class... As, class B0, class... Bs>
      struct is_compatible<Types<A0,As...>, Types<B0,Bs...>>
        : and_<is_compatible<A0, B0>, is_compatible<Types<As...>, Types<Bs...>>> {};
    }

    template <class A, class B>
    using IsCompatible = detail::is_compatible<A,B>;

  } // end namespace traits
  
  namespace concepts
  {
    /// Types support A + B
    template <class A, class B>
    using Addable = traits::IsAddable<A, B>;
    
    /// Types support A * B
    template <class A, class B>
    using Multiplicable = traits::IsMultiplicable<A, B>;
    
    /// Types are the same, up to decay of qualifiers
    template <class A, class B>
    using Compatible = traits::IsCompatible<A, B>;
    
  } // end namespace concepts
  
} // end namespace AMDiS
