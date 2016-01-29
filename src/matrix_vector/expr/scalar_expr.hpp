/** \file ScalarExpr.hpp */

#pragma once

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

#include "base_expr.hpp" // for BaseExpr

namespace AMDiS
{
  /// Expression to encapsulate scalars
  template <class Value>
  struct ScalarExpr
    : public BaseExpr<ScalarExpr<Value>>
  {
    using Self       = ScalarExpr;
    using expr_base  = BaseExpr<Self>;
    using size_type  = small_t;
    using value_type = Value;

    static constexpr int _SIZE = 1;
    static constexpr int _ROWS = 1;
    static constexpr int _COLS = 1;

  public:
    /// construcor takes the factor \p factor_
    constexpr ScalarExpr(value_type factor_)
      : factor(factor_) {}

    /// access the elements of an expr.
    constexpr value_type operator()(size_type = 0, size_type = 0) const
    {
      return factor;
    }

    /// cast operator for assignment
    operator value_type() const
    {
      return factor;
    }

  private:
    value_type factor;
  };

  /// Size of ScalarExpr
  template <class Value>
  inline size_t size(ScalarExpr<Value> const&)
  {
    return 1;
  }

  /// Size of ScalarExpr
  template <class Value>
  inline size_t num_rows(ScalarExpr<Value> const&)
  {
    return 1;
  }

  /// Size of ScalarExpr
  template <class Value>
  inline size_t num_cols(ScalarExpr<Value> const&)
  {
    return 1;
  }

  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    template <class V>
    struct category<ScalarExpr<V>>
    {
      using tag = typename category<V>::tag;
      using value_type = V;
      using size_type = small_t;
    };
    /// \endcond
  }

} // end namespace AMDiS
