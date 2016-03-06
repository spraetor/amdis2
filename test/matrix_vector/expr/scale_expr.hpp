/** \file ScaleExpr.hpp */

#pragma once

#include <operations/functors.hpp>

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>
#include <traits/scalar_types.hpp>
#include <traits/traits.hpp>

#include "matrix_vector/expr/shaped_expr.hpp"

namespace AMDiS
{
  /// Expression with one argument
  template <class Value, class Expr, bool from_left, class Functor>
  struct ScaleExpr
    : public ShapedExpr_t<Expr, ScaleExpr<Value, Expr, from_left, Functor>>
  {
    using Self       = ScaleExpr;

    using value_type = Value_t<Expr>;
    using size_type  = Size_t<Expr>;

    static constexpr int _SIZE = Expr::_SIZE;
    static constexpr int _ROWS = Expr::_ROWS;
    static constexpr int _COLS = Expr::_COLS;

  public:
    /// constructor takes the factor \p v and and expression \p A
    ScaleExpr(Value v, Expr const& A)
      : value(v), expr(&A)
    {
      MSG("ScaleExpr(" << v << ", " << &A << ")");
    }

    ~ScaleExpr() { MSG("~ScaleExpr()"); }

    ScaleExpr(Self const&) = delete;
//     ScaleExpr(Self&&)      = default;

    Self& operator=(Self const&) = delete;
//     Self& operator=(Self&&)      = default;

    /// move constructor
    ScaleExpr(Self&& other)
      : value(other.value),
        expr(other.expr)
    {
      other.expr = NULL;
    }

    Self& operator=(Self&& other)
    {
      delete expr;
      value = other.value;
      expr = other.expr;

      other.expr = NULL;
      return *this;
    }

    /// access the elements of an expr.
    value_type operator()(size_type i) const
    {
      return apply(i, bool_<from_left>{});
    }

    /// access the elements of a matrix-expr.
//     value_type operator()(size_type i, size_type j) const
//     {
//       return apply(i, j, bool_<from_left>());
//     }

    Expr const& get_first() const
    {
      return *expr;
    }

  private:
    // scale from left
    value_type apply(size_type i, true_) const
    {
      return fct( value, (*expr)(i) );
    }

    // scale from right
    value_type apply(size_type i, false_) const
    {
      return fct( (*expr)(i), value );
    }

    // scale from left
//     value_type apply(size_type i, size_type j, true_) const
//     {
//       return fct( value, expr(i,j) );
//     }
//
//     // scale from right
//     value_type apply(size_type i, size_type j, false_) const
//     {
//       return fct( expr(i,j), value );
//     }

  private:
    Value value;
    Expr const* expr;

    Functor fct;
  };


  /// Size of ScaleExpr
  template <class V, class E, bool l, class F>
  inline size_t size(ScaleExpr<V, E, l, F> const& expr)
  {
    size_t s = size(expr.get_first());
    std::cout << "size(ScaleExpr(" << &expr << ")) = " << s << "\n";
    return s;
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
      ScaleExpr<Value, E, true, functors::multiplies> >;

  // V * s
  template <class Value, class E>
  using RightScaleExpr
    = Requires_t<and_<concepts::Multiplicable<Value_t<E>, Value>,
                      concepts::Arithmetic<Value>>,
      ScaleExpr<Value, E, false, functors::multiplies> >;

  // V / s
  template <class Value, class E>
  using RightDivideExpr
    = Requires_t<and_<concepts::Multiplicable<Value_t<E>, Value>,
                      concepts::Arithmetic<Value>>,
      ScaleExpr<Value, E, false, functors::divides<Value_t<E>, Value>> >;

} // end namespace AMDiS
