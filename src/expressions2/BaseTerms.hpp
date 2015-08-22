/** \file BaseTerms.hpp */

#pragma once

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

namespace AMDiS 
{
  template <class Model>
  struct BaseTerm
  {
    Model& sub() { return static_cast<Model&>(*this); }
    Model const& sub() const { return static_cast<Model const&>(*this); }
  };
  
  template <class Model>
  struct VectorTerm : BaseTerm<Model> {};
  
  template <class Model>
  struct MatrixTerm : BaseTerm<Model> {};
  
  namespace traits 
  {
    /// \cond HIDDEN_SYMBOLS  
    // define categories for the terms
    template <class M> struct category<BaseTerm<M> > : category<M> {}; 
    template <class M> struct category<VectorTerm<M> > : category<M> {};
    template <class M> struct category<MatrixTerm<M> > : category<M> {};
    /// \endcond
  }
  
  // some elementary traits functions: size, num_rows, num_cols
  
  template <class T>
  inline size_t size(const VectorTerm<T>& t) { return size(t.sub()); }
  
  template <class T>
  inline size_t size(const MatrixTerm<T>& t) { return size(t.sub()); }
  
  template <class T>
  inline size_t num_rows(const VectorTerm<T>& t) { return num_rows(t.sub()); }
  
  template <class T>
  inline size_t num_rows(const MatrixTerm<T>& t) { return num_rows(t.sub()); }  
  
  template <class M>
  inline size_t num_cols(const VectorTerm<M>& t) { return num_cols(t.sub()); }
  
  template <class M>
  inline size_t num_cols(const MatrixTerm<M>& t) { return num_cols(t.sub()); }
  
  
  // determine shape of expression
  template <class Sub, class Model, class Enable = void>
  struct ShapedTerm 
  {
    using type = if_then_else< traits::is_vector<Sub>::value, 
                    VectorTerm<Model>,
                    if_then_else< traits::is_matrix<Sub>::value, 
                        MatrixTerm<Model>,
                        BaseTerm<Model> > >;
  };
  
  template <class Sub, class Model>
  using ShapedTerm_t = typename ShapedTerm<Sub, Model>::type;
  
} // end namespace AMDiS
