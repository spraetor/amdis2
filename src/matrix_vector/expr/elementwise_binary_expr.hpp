/** \file ElementwiseBinaryExpr.hpp */

#pragma once

#include <operations/functors.hpp>
#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>
#include <traits/traits.hpp>

#include "base_expr.hpp" // for ShapedExpr
#include <Math.h>

namespace AMDiS
{

  /// Expression with two arguments
  template <class E1, class E2, class Functor>
  struct ElementwiseBinaryExpr
    : public ShapedExpr_t<E1, ElementwiseBinaryExpr<E1, E2, Functor>>
  {
    using Self       = ElementwiseBinaryExpr;
    using expr_base  = ShapedExpr_t<E1, Self>;

    using value_type = typename std::result_of<Functor(Value_t<E1>, Value_t<E2>)>::type;
    using size_type  = traits::MaxSizeType<E1,E2>;
    using expr1_type = E1;
    using expr2_type = E2;

    constexpr static int _SIZE = math::max(E1::_SIZE, E2::_SIZE);
    constexpr static int _ROWS = math::max(E1::_ROWS, E2::_ROWS);
    constexpr static int _COLS = math::max(E1::_COLS, E2::_COLS);

    /// constructor takes two expressions
    ElementwiseBinaryExpr(expr1_type const& A, expr2_type const& B)
      : expr1(A), expr2(B)
    {
      TEST_EXIT_DBG( size(A) == size(B) )("Sizes do not match!\n");
    }

    /// access the elements of an expr.
    constexpr value_type operator()(size_type i) const
    {
      return fct( expr1(i), expr2(i) );
    }

    /// access the elements of a matrix-expr.
    constexpr value_type operator()(size_type i, size_type j) const
    {
      return fct( expr1(i, j), expr2(i, j) );
    }

    expr1_type const& get_first() const
    {
      return expr1;
    }
    expr2_type const& get_second() const
    {
      return expr2;
    }

  private:
    expr1_type const& expr1;
    expr2_type const& expr2;

    Functor fct;
  };


  /// Size of ElementwiseBinaryExpr
  template <class E1, class E2, class F>
  inline size_t size(ElementwiseBinaryExpr<E1, E2, F> const& expr)
  {
    return size(expr.get_first());
  }

  /// number of rows of ElementwiseBinaryExpr
  template <class E1, class E2, class F>
  inline size_t num_rows(ElementwiseBinaryExpr<E1, E2, F> const& expr)
  {
    return num_rows(expr.get_first());
  }

  /// number of columns of ElementwiseBinaryExpr
  template <class E1, class E2, class F>
  inline size_t num_cols(ElementwiseBinaryExpr<E1, E2, F> const& expr)
  {
    return num_cols(expr.get_first());
  }

  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    template <class E1, class E2, class F>
    struct category<ElementwiseBinaryExpr<E1,E2,F>>
    {
      using tag        = typename category<E1>::tag;
      using value_type = Value_t<ElementwiseBinaryExpr<E1,E2,F>>;
      using size_type  = MaxSizeType<E1,E2>;
    };
    /// \endcond
    
  } // end namespace traits

  // E1 + E2
  template <class E1, class E2>
  using PlusExpr 
    = ElementwiseBinaryExpr<E1, E2, functors::plus<Value_t<E1>, Value_t<E2>>>;

  // E1 - E2
  template <class E1, class E2>
  using MinusExpr 
    = ElementwiseBinaryExpr<E1, E2, functors::minus<Value_t<E1>, Value_t<E2>>>;


} // end namespace AMDiS
