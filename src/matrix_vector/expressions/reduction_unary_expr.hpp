/** \file ReductionUnaryExpr.hpp */

#pragma once


#include "traits/basic.hpp"
#include "base_expr.hpp" // for base_expr
#include "operations/reduction_functors.hpp"

namespace AMDiS {

  /// Expression with one argument, that reduces to a scalar
  template <class E, class Functor>
  struct ReductionUnaryExpr 
      : public BaseExpr< ReductionUnaryExpr<E, Functor> >
  {
    typedef ReductionUnaryExpr                   self;
    typedef BaseExpr<self>                  expr_base;
    
    typedef typename Functor::result_type  value_type;
    typedef typename E::size_type           size_type;
    typedef E                               expr_type;
    
    static constexpr int _SIZE = 1;
    static constexpr int _ROWS = 1;
    static constexpr int _COLS = 1;
    
  private:
    static constexpr int ARG_SIZE = E::_SIZE;
    
  public:
    /// constructor takes on expression \p A.
    ReductionUnaryExpr(expr_type const& A) 
      : expr(A)
    { }
    
    /// access the elements of an expr.
    inline value_type operator()(size_type = 0, size_type = 0) const
    {
      return reduce(int_<ARG_SIZE>()) ;
    }
    
    /// cast operator for assignment to scalar.
    operator value_type() const
    {
      return reduce(int_<ARG_SIZE>());
    }
    
    expr_type const& get_first() const { return expr; }
    
  protected:
    template <int N>
    value_type reduce(int_<N>) const
    {
      using meta::FOR;
      value_type erg; Functor::init(erg);
      
      for (size_type i = 0; i < size(expr); ++i)
	Functor::update(erg, expr(i));
      
//       FOR<0,N>::accumulate(expr, erg, Functor());
      return Functor::post_reduction(erg);
    }
    
    value_type reduce(int_<-1>) const
    {
      value_type erg; Functor::init(erg);
      for (size_type i = 0; i < size(expr); ++i)
	Functor::update(erg, expr(i));
      return Functor::post_reduction(erg);
    }
    
  private:
    expr_type const& expr;
  };
  
  /// Size of ReductionUnaryExpr
  template <class E, class F>
  size_t size(ReductionUnaryExpr<E,F> const&) { return 1; }
  
  /// Size of ReductionUnaryExpr
  template <class E, class F>
  size_t num_rows(ReductionUnaryExpr<E,F> const&) { return 1; }
  
  /// Size of ReductionUnaryExpr
  template <class E, class F>
  size_t num_cols(ReductionUnaryExpr<E,F> const&) { return 1; }
  
  namespace traits {
    
    /// \cond HIDDEN_SYMBOLS
    template <class E, class F>
    struct category<ReductionUnaryExpr<E,F> > 
    {
      typedef typename category<typename F::result_type>::tag    tag;
      typedef typename F::result_type                     value_type;
      typedef typename ReductionUnaryExpr<E,F>::size_type  size_type;
    };
    /// \endcond
  }
  
  
  // norm |V|_1
  template <class E>
  struct OneNormExpr {
    typedef ReductionUnaryExpr<E, functors::one_norm_functor<typename E::value_type> > type;
  };
    
  // norm |V|_2
  template <class E>
  struct TwoNormExpr {
    typedef ReductionUnaryExpr<E, functors::two_norm_functor<typename E::value_type> > type;
  };
    
  // V*V
  template <class E>
  struct UnaryDotExpr {
    typedef ReductionUnaryExpr<E, functors::unary_dot_functor<typename E::value_type> > type;
  };
    
  // max(V)
  template <class E>
  struct MaxExpr {
    typedef ReductionUnaryExpr<E, functors::max_reduction_functor<typename E::value_type> > type;
  };
    
  // abs_max(V)
  template <class E>
  struct AbsMaxExpr {
    typedef ReductionUnaryExpr<E, functors::abs_max_reduction_functor<typename E::value_type> > type;
  };
    
  // min(V)
  template <class E>
  struct MinExpr {
    typedef ReductionUnaryExpr<E, functors::min_reduction_functor<typename E::value_type> > type;
  };
    
  // max(V)
  template <class E>
  struct AbsMinExpr {
    typedef ReductionUnaryExpr<E, functors::abs_min_reduction_functor<typename E::value_type> > type;
  };
    
  // sum(V)
  template <class E>
  struct SumExpr {
    typedef ReductionUnaryExpr<E, functors::sum_reduction_functor<typename E::value_type> > type;
  };
    
  // prod(V)
  template <class E>
  struct ProdExpr {
    typedef ReductionUnaryExpr<E, functors::prod_reduction_functor<typename E::value_type> > type;
  };
  
} // end namespace AMDiS
