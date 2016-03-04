/** \file ScaleExpr.hpp */

#pragma once

#include <operations/functors.hpp>

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>
#include <traits/scalar_types.hpp>
#include <traits/traits.hpp>

#include "matrix_vector/expr/base_expr.hpp" // for shaped_expr

namespace AMDiS
{
  /// Expression with one argument
  template <class Value, class E, bool from_left, class Functor>
  struct ScaleExpr
    : public ShapedExpr_t<E, ScaleExpr<Value, E, from_left, Functor>>
  {
    using Self       = ScaleExpr;
    using expr_base  = ShapedExpr_t<E, Self>;

    using value_type = Value_t<E>;
    using size_type  = Size_t<E>;
    using expr_type  = E;

    static constexpr int _SIZE = E::_SIZE;
    static constexpr int _ROWS = E::_ROWS;
    static constexpr int _COLS = E::_COLS;

  public:
    /// constructor takes the factor \p v and and expression \p A
    constexpr ScaleExpr(Value v, expr_type const& A)
      : value(v), expr(A)
    {}

    /// access the elements of an expr.
    constexpr value_type operator()(size_type i) const
    {
      return apply(i, bool_<from_left>());
    }

    /// access the elements of a matrix-expr.
    constexpr value_type operator()(size_type i, size_type j) const
    {
      return apply(i, j, bool_<from_left>());
    }

    expr_type const& get_first() const
    {
      return expr;
    }

  private:
    // scale from left
    constexpr value_type apply(size_type i, true_) const
    {
      return fct( value, expr(i) );
    }

    // scale from right
    constexpr value_type apply(size_type i, false_) const
    {
      return fct( expr(i), value );
    }

    // scale from left
    constexpr value_type apply(size_type i, size_type j, true_) const
    {
      return fct( value, expr(i,j) );
    }

    // scale from right
    constexpr value_type apply(size_type i, size_type j, false_) const
    {
      return fct( expr(i,j), value );
    }

  private:
    Value value;
    expr_type const& expr;

    Functor fct;
  };


  /// Size of ScaleExpr
  template <class V, class E, bool l, class F>
  inline size_t size(ScaleExpr<V, E, l, F> const& expr)
  {
    return size(expr.get_first());
  }

  /// Size of ScaleExpr
  template <class V, class E, bool l, class F>
  inline size_t num_rows(ScaleExpr<V, E, l, F> const& expr)
  {
    return num_rows(expr.get_first());
  }

  /// Size of ScaleExpr
  template <class V, class E, bool l, class F>
  inline size_t num_cols(ScaleExpr<V, E, l, F> const& expr)
  {
    return num_cols(expr.get_first());
  }

  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    template <class V, class E, bool l, class F>
    struct category<ScaleExpr<V,E,l,F>> : category<E> {};
    /// \endcond

  } // end namespace traits


  // s * V
  template <class Value, class E>
  using LeftScaleExpr
    = Requires_t<and_<concepts::Multiplicable<Value_t<E>, Value>,
                      concepts::Arithmetic<Value>>,
      ScaleExpr<Value, E, true, functors::multiplies<Value_t<E>, Value>> >;

  // V * s
  template <class Value, class E>
  using RightScaleExpr
    = Requires_t<and_<concepts::Multiplicable<Value_t<E>, Value>,
                      concepts::Arithmetic<Value>>,
      ScaleExpr<Value, E, false, functors::multiplies<Value_t<E>, Value>> >;

  // V / s
  template <class Value, class E>
  using RightDivideExpr
    = Requires_t<and_<concepts::Multiplicable<Value_t<E>, Value>,
                      concepts::Arithmetic<Value>>,
      ScaleExpr<Value, E, false, functors::divides<Value_t<E>, Value>> >;

} // end namespace AMDiS
