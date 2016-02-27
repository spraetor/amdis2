#pragma once

// AMDiS includes
#include "traits/basic.hpp"
#include "traits/meta_basic.hpp"
#include "traits/traits_fwd.hpp"

namespace AMDiS
{
  template <class Model>
  struct BaseTerm
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
  struct VectorTerm : BaseTerm<Model> {};

  template <class Model>
  struct MatrixTerm : BaseTerm<Model> {};

  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    // define categories for the terms
    template <class M> struct category<BaseTerm<M>> : category<M> {};
    template <class M> struct category<VectorTerm<M>> : category<M> {};
    template <class M> struct category<MatrixTerm<M>> : category<M> {};
    /// \endcond
    
  } // end namespace traits

  // some elementary traits functions: size, num_rows, num_cols

  template <class T>
  inline size_t size(VectorTerm<T> const& t)
  {
    return size(t.sub());
  }

  template <class T>
  inline size_t size(MatrixTerm<T> const& t)
  {
    return size(t.sub());
  }

  template <class T>
  inline size_t num_rows(VectorTerm<T> const& t)
  {
    return num_rows(t.sub());
  }

  template <class T>
  inline size_t num_rows(MatrixTerm<T> const& t)
  {
    return num_rows(t.sub());
  }

  template <class M>
  inline size_t num_cols(VectorTerm<M> const& t)
  {
    return num_cols(t.sub());
  }

  template <class M>
  inline size_t num_cols(MatrixTerm<M> const& t)
  {
    return num_cols(t.sub());
  }


  // determine shape of expression
  template <class Sub, class Model>
  struct ShapedTerm
  {
    using type 
      = if_then_else< traits::is_vector<Sub>::value,    VectorTerm<Model>,
        if_then_else< traits::is_matrix<Sub>::value,    MatrixTerm<Model>,
                                                        BaseTerm<Model> > >;
  };

  template <class Sub, class Model>
  using ShapedTerm_t = typename ShapedTerm<Sub, Model>::type;

} // end namespace AMDiS
