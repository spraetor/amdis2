/** \file BaseExpr.hpp */

#pragma once

#include <matrix_vector/Forward.h>
#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

namespace AMDiS
{
  template <class Model>
  struct BaseExpr
  {
    Model& sub()
    {
      return static_cast<Model&>(*this);
    }
    Model const& sub() const
    {
      return static_cast<Model const&>(*this);
    }
  };

  template <class Model>
  struct VectorExpr : BaseExpr<Model> {};

  template <class Model>
  struct MatrixExpr : BaseExpr<Model> {};

  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    // define categories for the expressions
    template <class M> struct category<BaseExpr<M>> : category<M> {};
    template <class M> struct category<VectorExpr<M>> : category<M> {};
    template <class M> struct category<MatrixExpr<M>> : category<M> {};
    /// \endcond
  }


  template <class M>
  inline size_t num_cols(const VectorExpr<M>& t)
  {
    return num_cols(t.sub());
  }
  template <class M>
  inline size_t num_cols(const MatrixExpr<M>& t)
  {
    return num_cols(t.sub());
  }

  template <class T>
  inline size_t num_rows(const VectorExpr<T>& t)
  {
    return num_rows(t.sub());
  }
  template <class T>
  inline size_t num_rows(const MatrixExpr<T>& t)
  {
    return num_rows(t.sub());
  }

  template <class T>
  inline size_t size(const VectorExpr<T>& t)
  {
    return size(t.sub());
  }
  template <class T>
  inline size_t size(const MatrixExpr<T>& t)
  {
    return size(t.sub());
  }


  // determine shape of expression
  template <class Sub, class Model, class Enable = void>
  struct ShapedExpr
  {
    using type = if_then_else<
                 traits::is_vector<Sub>::value,
                 VectorExpr<Model>,
                 if_then_else<traits::is_matrix<Sub>::value,
                 MatrixExpr<Model>,
                 BaseExpr<Model>> >;
  };

  template <class Sub, class Model>
  using ShapedExpr_t = typename ShapedExpr<Sub, Model>::type;

} // end namespace AMDiS
