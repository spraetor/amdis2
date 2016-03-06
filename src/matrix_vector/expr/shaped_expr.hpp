#pragma once

#include "traits/basic.hpp"
#include "matrix_vector/expr/base_expr.hpp"
#include "matrix_vector/ExprConcepts.hpp"

namespace AMDiS
{
  // TODO: write begin(expr) and end(expr), i.e. provide general
  // iterator for expressions

  // determine shape of expression
  template <class Sub, class Model>
  struct ShapedExpr
  {
    using type
      = if_then_else< concepts::VectorExpression<Sub>::value,    VectorExpr<Model>,
        if_then_else< concepts::MatrixExpression<Sub>::value,    MatrixExpr<Model>,
                                                                 BaseExpr<Model> > >;
  };

  template <class Sub, class Model>
  using ShapedExpr_t = typename ShapedExpr<Sub, Model>::type;

} // end namespace AMDiS
