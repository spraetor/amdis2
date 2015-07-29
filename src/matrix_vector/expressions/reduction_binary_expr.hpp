/** \file ReductionBinaryExpr.hpp */

#pragma once

#include "traits/basic.hpp"
#include "traits/category.hpp"
#include "base_expr.hpp" // for base_expr

namespace AMDiS {

  /// Expression with two arguments, that reduces to a scalar
  template <class E1, class E2, class Functor>
  struct ReductionBinaryExpr 
      : public BaseExpr< ReductionBinaryExpr<E1, E2, Functor> >
  {
    typedef ReductionBinaryExpr                               self;
    typedef BaseExpr<self>                               expr_base;
    
    typedef typename Functor::result_type               value_type;
    typedef typename traits::max_size_type<E1,E2>::type  size_type;
    typedef E1                                          expr1_type;
    typedef E2                                          expr2_type;
    
    static constexpr int _SIZE = 1;
    static constexpr int _ROWS = 1;
    static constexpr int _COLS = 1;
    
  private:
    static constexpr int ARG_SIZE = MAX(E1::_SIZE, E2::_SIZE);
    
  public:
    /// constructor takes two expression \p A and \p B.
    ReductionBinaryExpr(expr1_type const& A, expr2_type const& B) 
	: expr1(A), expr2(B)
    { 
      TEST_EXIT_DBG( size(A) == size(B) )("Sizes do not match!\n");
    }
    
    /// access the elements of an expr.
    value_type operator()(size_type = 0, size_type = 0) const
    {
      return reduce(int_<ARG_SIZE>());
    }
    
    /// cast operator for assignment to scalar.
    operator value_type() const
    {
      return reduce(int_<ARG_SIZE>());
    }
    
    expr1_type const& get_first() const { return expr1; }
    expr2_type const& get_second() const { return expr2; }
    
  protected:
    template <int N>
    value_type reduce(int_<N>) const
    {
      using meta::FOR;
      value_type erg; Functor::init(erg);
      for (size_type i = 0; i < N; ++i)
	Functor::update(erg, expr1(i), expr2(i));
//       FOR<0,N>::inner_product(expr1, expr2, erg, Functor());
      return Functor::post_reduction(erg);
    }
    
    value_type reduce(int_<-1>) const
    {
      value_type erg; Functor::init(erg);
      for (size_type i = 0; i < size(expr1); ++i)
	Functor::update(erg, expr1(i), expr2(i));
      return Functor::post_reduction(erg);
    }
    
  private:
    expr1_type const& expr1;
    expr2_type const& expr2;
  };
  
  /// Size of ReductionBinaryExpr
  template <class E1, class E2, class F>
  size_t size(ReductionBinaryExpr<E1,E2,F> const&) { return 1; }
  
  /// Size of ReductionBinaryExpr
  template <class E1, class E2, class F>
  size_t num_rows(ReductionBinaryExpr<E1,E2,F> const&) { return 1; }
  
  /// Size of ReductionBinaryExpr
  template <class E1, class E2, class F>
  size_t num_cols(ReductionBinaryExpr<E1,E2,F> const&) { return 1; }
  
  namespace traits {
    
    /// \cond HIDDEN_SYMBOLS
    template <class E1, class E2, class F>
    struct category<ReductionBinaryExpr<E1,E2,F> > 
    {
      typedef typename category<typename F::result_type>::tag         tag;
      typedef typename F::result_type                          value_type;
      typedef typename ReductionBinaryExpr<E1,E2,F>::size_type  size_type;
    };
    /// \endcond
  }
  
  // standard inner product
  template <class E1, class E2>
  struct DotExpr {
    typedef ReductionBinaryExpr<E1, E2, 
	      functors::dot_functor<typename E1::value_type, typename E2::value_type> > type;
  };
  
} // end namespace AMDiS
