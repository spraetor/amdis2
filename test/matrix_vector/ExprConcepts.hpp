#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// AMDiS headers
#include "matrix_vector/expr/base_expr.hpp"
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
    struct Expr<T> : std::is_base_of< BaseExpr<Decay_t<T>>, Decay_t<T> > {};
//     {
//     private:
//       using Type = Decay_t<T>;
//
//       template <class T1> static Ints<T1::_SIZE, T1::_ROWS, T1::_COLS> test1(int);
//       template <class>    static void test1(...);
//
//       static constexpr bool value0 = traits::HasValueType<Type>::value &&
//                                      traits::HasSizeType<Type>::value;
//       static constexpr bool value1 = !std::is_void<decltype(test1<Type>(0))>::value;
//       //static constexpr bool value2 = std::is_base_of<BaseExpr<Type>, Type>::value;
//
//     public:
//       static constexpr bool value = value0 && value1; // && value2;
//     };


    /// Expressions concept for vectors
    template <class... Ts>
    struct VectorExpression : and_<VectorExpression<Ts>...> {};

    template <class T>
    struct VectorExpression<T> : std::is_base_of< VectorExpr<Decay_t<T>>, Decay_t<T> > {};


    /// Expressions concept for matrices
    template <class... Ts>
    struct MatrixExpression : and_<VectorExpression<Ts>...> {};

    template <class T>
    struct MatrixExpression<T> : std::is_base_of< MatrixExpr<Decay_t<T>>, Decay_t<T> > {};


  } // end namespace concepts


  namespace requires
  {
    /// Test whether one of the arguments is an expression
    template <class Result, class... Args>
    using Expr = Requires_t<or_<concepts::Expr<Args>...>, Result>;

  } // end namespace require
} // end namespace AMDiS
