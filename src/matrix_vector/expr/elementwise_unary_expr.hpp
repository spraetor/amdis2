/** \file ElementwiseUnaryExpr.hpp */

#pragma once

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

#include "base_expr.hpp" // for ShapedExpr

namespace AMDiS {

  /// Expression with one argument
  template <class E, class Functor>
  struct ElementwiseUnaryExpr 
      : public ShapedExpr<E, ElementwiseUnaryExpr<E, Functor> >::type
  {
    typedef ElementwiseUnaryExpr                     Self;
    typedef typename ShapedExpr<E, Self>::type  expr_base;
    
    typedef Result_t<Functor>                  value_type;
    typedef Size_t<E>                           size_type;
    typedef E                                   expr_type;
    
    static constexpr int _SIZE = E::_SIZE;
    static constexpr int _ROWS = E::_ROWS;
    static constexpr int _COLS = E::_COLS;
    
    /// constructor takes an expression
    ElementwiseUnaryExpr(expr_type const& A) 
	: expr(A) 
    { }
    
    /// access the elements of an expr.
    value_type operator()(size_type i) const
    { 
      return Functor::apply( expr(i) );
    }
    
    /// access the elements of a matrix-expr.
    value_type operator()(size_type i, size_type j) const
    { 
      return Functor::apply( expr(i, j) );
    }
    
    expr_type const& get_first() const { return expr; }
    
  private:
    expr_type const& expr;
  };
  
  
  /// Size of ElementwiseUnaryExpr
  template <class E, class F>
  size_t size(ElementwiseUnaryExpr<E, F> const& expr)
  {
    return size(expr.get_first());
  }
  
  /// number of rows of ElementwiseUnaryExpr
  template <class E, class F>
  size_t num_rows(ElementwiseUnaryExpr<E, F> const& expr)
  {
    return num_rows(expr.get_first());
  }
  
  /// number of columns of ElementwiseUnaryExpr
  template <class E, class F>
  size_t num_cols(ElementwiseUnaryExpr<E, F> const& expr)
  {
    return num_cols(expr.get_first());
  }
  
  namespace traits {
    
    /// \cond HIDDEN_SYMBOLS
    template <class M, class F>
    struct category<ElementwiseUnaryExpr<M,F> > : category<M> {};
    /// \endcond
  }
} // end namespace AMDiS
