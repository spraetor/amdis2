/** \file ScaleExpr.hpp */

#pragma once

#include <boost/numeric/mtl/operation/sfunctor.hpp>

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

#include "base_expr.hpp" // for shaped_expr

namespace AMDiS {

  /// Expression with one argument
  template <class Value, class E, bool from_left, class Functor>
  struct ScaleExpr 
      : public ShapedExpr<E, ScaleExpr<Value, E, from_left, Functor> >::type
  {
    typedef ScaleExpr                                Self;
    typedef typename ShapedExpr<E, Self>::type  expr_base;
    
    typedef Value_t<E>                         value_type;
    typedef Size_t<E>                           size_type;
    typedef E                                   expr_type;
    
    static constexpr int _SIZE = E::_SIZE;
    static constexpr int _ROWS = E::_ROWS;
    static constexpr int _COLS = E::_COLS;
    
  public:
    /// constructor takes the factor \p v and and expression \p A
    ScaleExpr(Value v, expr_type const& A) 
	: value(v), expr(A) 
    { }
    
    /// access the elements of an expr.
    value_type operator()(size_type i) const
    { 
      return apply(i, bool_<from_left>());
    }
    
    /// access the elements of a matrix-expr.
    value_type operator()(size_type i, size_type j) const
    { 
      return apply(i, j, bool_<from_left>());
    }
    
    expr_type const& get_first() const { return expr; }
    
  private:
    // scale from left
    value_type apply(size_type i, true_) const 
    {
      return Functor::apply( value, expr(i) );
    }
    
    // scale from right
    value_type apply(size_type i, false_) const 
    {
      return Functor::apply( expr(i), value );
    }
    
    // scale from left
    value_type apply(size_type i, size_type j, true_) const 
    {
      return Functor::apply( value, expr(i,j) );
    }
    
    // scale from right
    value_type apply(size_type i, size_type j, false_) const 
    {
      return Functor::apply( expr(i,j), value );
    }
    
  private:
    Value value;
    expr_type const& expr;
  };
  
  
  /// Size of ScaleExpr
  template <class V, class E, bool l, class F>
  size_t size(ScaleExpr<V, E, l, F> const& expr)
  {
    return size(expr.get_first());
  }
  
  /// Size of ScaleExpr
  template <class V, class E, bool l, class F>
  size_t num_rows(ScaleExpr<V, E, l, F> const& expr)
  {
    return num_rows(expr.get_first());
  }
  
  /// Size of ScaleExpr
  template <class V, class E, bool l, class F>
  size_t num_cols(ScaleExpr<V, E, l, F> const& expr)
  {
    return num_cols(expr.get_first());
  }
  
  namespace traits {
    
    /// \cond HIDDEN_SYMBOLS
    template <class V, class E, bool l, class F>
    struct category<ScaleExpr<V,E,l,F> > : category<E> {};
    /// \endcond
  }

  // s * V
  template <class Value, class E>
  struct LeftScaleExpr : enable_if< 
    and_< traits::is_multiplicable<typename E::value_type, Value>,
	  traits::is_scalar<Value> >,
    ScaleExpr<Value, E, true, mtl::sfunctor::times<typename E::value_type, Value> > > {};
    
  // V * s
  template <class Value, class E>
  struct RightScaleExpr : enable_if< 
    and_< traits::is_multiplicable<typename E::value_type, Value>,
	  traits::is_scalar<Value> >,
    ScaleExpr<Value, E, false, mtl::sfunctor::times<typename E::value_type, Value> > > {};
  
  // V / s
  template <class Value, class E>
  struct RightDivideExpr : enable_if< 
    and_< traits::is_multiplicable<typename E::value_type, Value>,
	  traits::is_scalar<Value> >,
    ScaleExpr<Value, E, false, mtl::sfunctor::divide<typename E::value_type, Value> > > {};
  
} // end namespace AMDiS
