/** \file BaseExpr.hpp */

#pragma once

#include "matrix_vector/Forward.hpp"
#include "traits/basic.hpp"
#include "traits/traits_fwd.hpp"

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
    template <class M> struct category<BaseExpr<M>>   : category<M> {};
    template <class M> struct category<VectorExpr<M>> : category<M> {};
    template <class M> struct category<MatrixExpr<M>> : category<M> {};
    /// \endcond

  } // end namespace traits


  template <class M>
  inline size_t num_cols(VectorExpr<M> const& t)
  {
    return num_cols(t.sub());
  }
  template <class M>
  inline size_t num_cols(MatrixExpr<M> const& t)
  {
    return num_cols(t.sub());
  }

  template <class T>
  inline size_t num_rows(VectorExpr<T> const& t)
  {
    return num_rows(t.sub());
  }
  template <class T>
  inline size_t num_rows(MatrixExpr<T> const& t)
  {
    return num_rows(t.sub());
  }

  template <class T>
  inline size_t size(VectorExpr<T> const& t)
  {
    size_t s = size(t.sub());
    std::cout << "size(VectorExpr(" << &t.sub() << " => " << typeid(t.sub()).name() << ") = " << s << "\n";
    return s;
  }

  template <class T>
  inline size_t size(MatrixExpr<T> const& t)
  {
    return size(t.sub());
  }

} // end namespace AMDiS
