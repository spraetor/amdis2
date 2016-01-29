/** \file ReductionUnaryExpr.hpp */

#pragma once


#include <traits/basic.hpp>
#include "base_expr.hpp" // for base_expr
#include <operations/reduction_functors.hpp>

namespace AMDiS
{
  /// Expression with one argument, that reduces to a scalar
  template <class E, class Functor>
  struct ReductionUnaryExpr
    : public BaseExpr<ReductionUnaryExpr<E, Functor>>
  {
    using Self       = ReductionUnaryExpr;
    using expr_base  = BaseExpr<Self>;

    using value_type = Result_t<Functor>;
    using size_type  = Size_t<E>;
    using expr_type  = E;

    static constexpr int _SIZE = 1;
    static constexpr int _ROWS = 1;
    static constexpr int _COLS = 1;

  private:
    static constexpr int ARG_SIZE = E::_SIZE;

  public:
    /// constructor takes on expression \p A.
    constexpr ReductionUnaryExpr(expr_type const& A)
      : expr(A)
    { }

    /// access the elements of an expr.
    value_type operator()(size_type = 0, size_type = 0) const
    {
      return reduce(int_<ARG_SIZE>()) ;
    }

    /// cast operator for assignment to scalar.
    operator value_type() const
    {
      return reduce(int_<ARG_SIZE>());
    }

    expr_type const& get_first() const
    {
      return expr;
    }

  protected:
    // reduce the expression, if length is known at compiletime
    template <int N>
    value_type reduce(int_<N>) const
    {
      using meta::FOR;
      value_type erg;
      Functor::init(erg);

      for (size_type i = 0; i < size(expr); ++i)
        Functor::update(erg, expr(i));

      //       FOR<0,N>::accumulate(expr, erg, Functor());
      return Functor::post_reduction(erg);
    }

    // reduce the expression, if length is known only at runtime
    value_type reduce(int_<-1>) const
    {
      value_type erg;
      Functor::init(erg);
      for (size_type i = 0; i < size(expr); ++i)
        Functor::update(erg, expr(i));
      return Functor::post_reduction(erg);
    }

  private:
    expr_type const& expr;
  };

  /// Size of ReductionUnaryExpr
  template <class E, class F>
  size_t size(ReductionUnaryExpr<E,F> const&)
  {
    return 1;
  }

  /// Size of ReductionUnaryExpr
  template <class E, class F>
  size_t num_rows(ReductionUnaryExpr<E,F> const&)
  {
    return 1;
  }

  /// Size of ReductionUnaryExpr
  template <class E, class F>
  size_t num_cols(ReductionUnaryExpr<E,F> const&)
  {
    return 1;
  }

  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    template <class E, class F>
    struct category<ReductionUnaryExpr<E,F>>
    {
      using tag = typename category<Result_t<F>>::tag;
      using value_type = Result_t<F>;
      using size_type  = Size_t<ReductionUnaryExpr<E,F>>;
    };
    /// \endcond
    
  } // end namespace traits


  // norm |V|_1
  template <class E>
  struct OneNormExpr
  {
    using type = ReductionUnaryExpr<E, functors::one_norm_functor<Value_t<E>>>;
  };

  // norm |V|_2
  template <class E>
  struct TwoNormExpr
  {
    using type = ReductionUnaryExpr<E, functors::two_norm_functor<Value_t<E>>>;
  };

  // V*V
  template <class E>
  struct UnaryDotExpr
  {
    using type = ReductionUnaryExpr<E, functors::unary_dot_functor<Value_t<E>>>;
  };

  // max(V)
  template <class E>
  struct MaxExpr
  {
    using type = ReductionUnaryExpr<E, functors::max_reduction_functor<Value_t<E>>>;
  };

  // abs_max(V)
  template <class E>
  struct AbsMaxExpr
  {
    using type = ReductionUnaryExpr<E, functors::abs_max_reduction_functor<Value_t<E>>>;
  };

  // min(V)
  template <class E>
  struct MinExpr
  {
    using type = ReductionUnaryExpr<E, functors::min_reduction_functor<Value_t<E>>>;
  };

  // max(V)
  template <class E>
  struct AbsMinExpr
  {
    using type = ReductionUnaryExpr<E, functors::abs_min_reduction_functor<Value_t<E>>>;
  };

  // sum(V)
  template <class E>
  struct SumExpr
  {
    using type = ReductionUnaryExpr<E, functors::sum_reduction_functor<Value_t<E>>>;
  };

  // prod(V)
  template <class E>
  struct ProdExpr
  {
    using type = ReductionUnaryExpr<E, functors::prod_reduction_functor<Value_t<E>>>;
  };

} // end namespace AMDiS
