/** \file ElementwiseUnaryExpr.hpp */

#pragma once

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

#include "base_expr.hpp" // for ShapedExpr

namespace AMDiS {

  /// Expression with one argument
  template <class E, class Functor>
  struct ElementwiseUnaryExpr 
      : public ShapedExpr_t<E, ElementwiseUnaryExpr<E, Functor> >
  {
    using Self = ElementwiseUnaryExpr;
    using expr_base = ShapedExpr_t<E, Self>;
    
    using value_type = typename std::result_of<Functor(Value_t<E>)>::type;
    using size_type  = Size_t<E>;
    using expr_type  = E;
    
    constexpr static int _SIZE = E::_SIZE;
    constexpr static int _ROWS = E::_ROWS;
    constexpr static int _COLS = E::_COLS;
    
    /// constructor takes an expression
    constexpr ElementwiseUnaryExpr(expr_type const& A) 
      : expr(A) 
    {}
    
    /// access the elements of an expr.
    constexpr value_type operator()(size_type i) const
    { 
      return fct( expr(i) );
    }
    
    /// access the elements of a matrix-expr.
    constexpr value_type operator()(size_type i, size_type j) const
    { 
      return fct( expr(i, j) );
    }
    
    expr_type const& get_first() const { return expr; }
    
  private:
    expr_type const& expr;
    
    Functor fct;
  };
  
  
  /// Size of ElementwiseUnaryExpr
  template <class E, class F>
  inline size_t size(ElementwiseUnaryExpr<E, F> const& expr)
  {
    return size(expr.get_first());
  }
  
  /// number of rows of ElementwiseUnaryExpr
  template <class E, class F>
  inline size_t num_rows(ElementwiseUnaryExpr<E, F> const& expr)
  {
    return num_rows(expr.get_first());
  }
  
  /// number of columns of ElementwiseUnaryExpr
  template <class E, class F>
  inline size_t num_cols(ElementwiseUnaryExpr<E, F> const& expr)
  {
    return num_cols(expr.get_first());
  }
  
  namespace traits 
  {
    /// \cond HIDDEN_SYMBOLS
    template <class M, class F>
    struct category<ElementwiseUnaryExpr<M,F> > : category<M> {};
    /// \endcond
  }
} // end namespace AMDiS
