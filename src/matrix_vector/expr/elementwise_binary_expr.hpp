/** \file ElementwiseBinaryExpr.hpp */

#pragma once

#include <boost/numeric/mtl/operation/sfunctor.hpp>

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

#include "base_expr.hpp" // for ShapedExpr
#include <Math.h>

namespace AMDiS {

  /// Expression with two arguments
  template <class E1, class E2, class Functor>
  struct ElementwiseBinaryExpr
      : public ShapedExpr<E1, ElementwiseBinaryExpr<E1, E2, Functor> >::type
  {
    typedef ElementwiseBinaryExpr                             Self;
    typedef typename ShapedExpr<E1, Self>::type          expr_base;
    
    typedef Result_t<Functor>                           value_type;
    typedef traits::max_size_type<E1,E2>                 size_type;
    typedef E1                                          expr1_type;
    typedef E2                                          expr2_type;
    
    static constexpr int _SIZE = math::max(E1::_SIZE, E2::_SIZE);
    static constexpr int _ROWS = math::max(E1::_ROWS, E2::_ROWS);
    static constexpr int _COLS = math::max(E1::_COLS, E2::_COLS);
    
    /// constructor takes two expressions
    ElementwiseBinaryExpr(expr1_type const& A, expr2_type const& B) 
	: expr1(A), expr2(B) 
    { 
      TEST_EXIT_DBG( size(A) == size(B) )("Sizes do not match!\n");
    }
    
    /// access the elements of an expr.
    value_type operator()(size_type i) const
    { 
      return Functor::apply( expr1(i), expr2(i) );
    }
    
    /// access the elements of a matrix-expr.
    value_type operator()(size_type i, size_type j) const
    { 
      return Functor::apply( expr1(i, j), expr2(i, j) );
    }
    
    expr1_type const& get_first() const { return expr1; }
    expr2_type const& get_second() const { return expr2; }
    
  private:
    expr1_type const& expr1;
    expr2_type const& expr2;
  };
  
  
  /// Size of ElementwiseBinaryExpr
  template <class E1, class E2, class F>
  size_t size(ElementwiseBinaryExpr<E1, E2, F> const& expr)
  {
    return size(expr.get_first());
  }
  
  /// number of rows of ElementwiseBinaryExpr
  template <class E1, class E2, class F>
  size_t num_rows(ElementwiseBinaryExpr<E1, E2, F> const& expr)
  {
    return num_rows(expr.get_first());
  }
  
  /// number of columns of ElementwiseBinaryExpr
  template <class E1, class E2, class F>
  size_t num_cols(ElementwiseBinaryExpr<E1, E2, F> const& expr)
  {
    return num_cols(expr.get_first());
  }
  
  namespace traits {
    
    /// \cond HIDDEN_SYMBOLS
    template <class E1, class E2, class F>
    struct category< ElementwiseBinaryExpr<E1,E2,F> > 
    {
      typedef typename category<E1>::tag          tag;
      typedef Result_t<F>                  value_type;
      typedef max_size_type<E1,E2>          size_type;
    };
    /// \endcond
  }

  // E1 + E2
  template <class E1, class E2>
  struct PlusExpr {
    typedef ElementwiseBinaryExpr<E1, E2, 
      mtl::sfunctor::plus<Value_t<E1>, Value_t<E2> > > type;
  };
      
  // E1 - E2
  template <class E1, class E2>
  struct MinusExpr {
    typedef ElementwiseBinaryExpr<E1, E2, 
      mtl::sfunctor::minus<Value_t<E1>, Value_t<E2> > > type;
  };
      
  
} // end namespace AMDiS
