#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// AMDiS headers
#include "traits/basic.hpp"
#include "traits/meta_basic.hpp"
#include "traits/concepts_base.hpp"

namespace AMDiS
{
  namespace concepts
  {
    /// Expressions concepts
    template <class... Ts>
    struct Expr : and_<Expr<Ts>...> {};

    template <class T>
    struct Expr<T>
    {
    private:
      using Type = Decay_t<T>;

      template <class T1> static Ints<T1::_SIZE, T1::_ROWS, T1::_COLS> test1(int);
      template <class>    static void test1(...);

      static constexpr bool value0 = traits::HasValueType<Type>::value && 
				     traits::HasSizeType<Type>::value;
      static constexpr bool value1 = not std::is_void<decltype(test1<Type>(0))>::value;
      static constexpr bool value2 = std::is_base_of<BaseExpr<Type>, Type>::value;

    public:
      static constexpr bool value = value0 && value1 && value2;
    };

  } // end namespace concepts


  namespace requires
  {
    /// Test whether one of the arguments is an expression
    template <class Result, class... Args>
    using Expr = Requires_t<or_<concepts::Expr<Args>...>, Result>;

  } // end namespace require
} // end namespace AMDiS
