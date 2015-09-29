/** \file ExprConcepts.hpp */

#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// AMDiS headers
#include <traits/basic.hpp>
#include <traits/meta_basic.hpp>
#include <traits/concepts_base.hpp>

namespace AMDiS
{
  namespace concepts
  {
    // expressions concepts
    template <class... Ts>
    struct Expr
    {
      constexpr static bool value = and_<Expr<Ts>...>::value;
    };

    template <class T>
    struct Expr<T>
    {
    private:
      using Type = typename std::decay<T>::type;

      template <class T1> static Ints<T1::_SIZE, T1::_ROWS, T1::_COLS> test1(int);
      template <class> static void test1(...);

      constexpr static bool value0 = traits::HasValueType<Type>::value && traits::HasSizeType<Type>::value;
      constexpr static bool value1 = !std::is_void<decltype(test1<Type>(0))>::value;
      constexpr static bool value2 = std::is_base_of<BaseExpr<Type>, Type>::value;

    public:
      static constexpr bool value = value0 && value1 && value2;
    };

  } // end namespace concepts


  namespace requires
  {
    // test whether one of the arguments is an expression
    template <class Result, class... Args>
    using Expr = Requires_t<or_<concepts::Expr<Args>...>, Result>;

  } // end namespace require
} // end namespace AMDiS
