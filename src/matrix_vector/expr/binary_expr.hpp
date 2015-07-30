/** \file binary_expr.hpp */

#pragma once


#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

#include "base_expr.hpp" // for shaped_expr

#include <operations/functors.hpp>
#include <traits/mult_type.hpp>
#include <Math.h>

namespace AMDiS {

  /// Expression with two arguments
  template <class E1, class E2, class Functor>
  struct VectorBinaryExpr
      : public VectorExpr< VectorBinaryExpr<E1, E2, Functor> >
  {
    typedef VectorBinaryExpr                                  Self;
    typedef VectorExpr<Self>                             expr_base;
    
    typedef Value_t<E1>                                 value_type;
    typedef traits::max_size_type<E1,E2>                 size_type;
    typedef E1                                          expr1_type;
    typedef E2                                          expr2_type;
    
    static constexpr int _SIZE = math::max(E1::_SIZE, E2::_SIZE);
    static constexpr int _ROWS = math::max(E1::_ROWS, E2::_ROWS);
    static constexpr int _COLS = math::max(E1::_COLS, E2::_COLS);
    
  public:
    /// constructor takes two expressions
    VectorBinaryExpr(expr1_type const& A, expr2_type const& B) 
	: expr1(A), expr2(B) 
    { 
      assert( size(A) == size(B) );
    }
    
    /// access the elements of an expr.
    value_type operator()(size_type i) const
    { 
      return Functor::apply( i, expr1, expr2 );
    }
    
    /// access the elements of a matrix-expr.
    value_type operator()(size_type i, size_type j) const
    { 
      return Functor::apply( i, j, expr1, expr2 );
    }
    
    expr1_type const& get_first() const { return expr1; }
    expr2_type const& get_second() const { return expr2; }
    
  private:
    expr1_type const& expr1;
    expr2_type const& expr2;
  };
  
  /// Size of VectorBinaryExpr
  template <class E1, class E2, class F>
  size_t size(VectorBinaryExpr<E1, E2, F> const& expr)
  {
    return size(expr.get_first());
  }
  
  /// number of rows of VectorBinaryExpr
  template <class E1, class E2, class F>
  size_t num_rows(VectorBinaryExpr<E1, E2, F> const& expr)
  {
    return num_rows(expr.get_first());
  }
  
  /// number of columns of VectorBinaryExpr
  template <class E1, class E2, class F>
  size_t num_cols(VectorBinaryExpr<E1, E2, F> const& expr)
  {
    return num_cols(expr.get_first());
  }
  
  namespace traits {
    
    /// \cond HIDDEN_SYMBOLS
    template <class E1, class E2, class F>
    struct category<VectorBinaryExpr<E1,E2,F> > 
    {
      typedef typename category<typename VectorBinaryExpr<E1,E2,F>::expr_base>::tag tag;
      typedef typename VectorBinaryExpr<E1,E2,F>::value_type  value_type;
      typedef typename VectorBinaryExpr<E1,E2,F>::size_type   size_type;
    };
    /// \endcond
  }
  
} // end namespace AMDiS
