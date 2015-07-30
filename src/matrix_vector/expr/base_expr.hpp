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
    Model& sub() { return static_cast<Model&>(*this); }
    Model const& sub() const { return static_cast<Model const&>(*this); }
  };
  
  template <class Model>
  struct VectorExpr : BaseExpr<Model> {};
  
  template <class Model>
  struct MatrixExpr : BaseExpr<Model> {};
  
  namespace traits 
  {
    /// \cond HIDDEN_SYMBOLS  
    // define categories for the expressions
    template <class M> struct category<BaseExpr<M> > : category<M> {}; 
    template <class M> struct category<VectorExpr<M> > : category<M> {};
    template <class M> struct category<MatrixExpr<M> > : category<M> {};
    
    // define element access for the expressions
    template <class M> struct at<BaseExpr<M> > : at<M> {};
    template <class M> struct at<VectorExpr<M> > : at<M> {};
    template <class M> struct at<MatrixExpr<M> > : at<M> {};
    
    // define size of an expression
    template <class M> struct size<BaseExpr<M> > : size<M> {};
    template <class M> struct size<VectorExpr<M> > : size<M> {};
    template <class M> struct size<MatrixExpr<M> > : size<M> {};
    
    // define number of rows of an expressions
    template <class M> struct num_rows<BaseExpr<M> > : num_rows<M> {};
    template <class M> struct num_rows<VectorExpr<M> > : num_rows<M> {};
    template <class M> struct num_rows<MatrixExpr<M> > : num_rows<M> {};
    
    // define number of columns for an expression
    template <class M> struct num_cols<BaseExpr<M> > : num_cols<M> {};
    template <class M> struct num_cols<VectorExpr<M> > : num_cols<M> {};
    template <class M> struct num_cols<MatrixExpr<M> > : num_cols<M> {};
    /// \endcond
  }
  
  
  
  template <class M>
  size_t inline num_cols(const VectorExpr<M>& t)
  {
      return num_cols(t.sub());
  }
  
  template <class T>
  size_t inline num_rows(const VectorExpr<T>& t)
  {
      return num_rows(t.sub());
  }
  
  template <class T>
  size_t inline size(const VectorExpr<T>& t)
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
	   if_then_else< traits::is_matrix<Sub>::value, 
			 MatrixExpr<Model>,
			 BaseExpr<Model> > >;
  };

  
} // end namespace AMDiS
