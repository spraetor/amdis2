/** \file BaseExpr.hpp */

#pragma once

#include "traits/basic.hpp"
#include "traits/category.hpp"
#include "traits/size.hpp"
#include "traits/num_rows.hpp"
#include "traits/num_cols.hpp"
#include "traits/at.hpp"

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
  
  namespace traits {
    
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
  
  // determine shape of expression
  template <class Sub, class Model, class Enable = void>
  struct ShapedExpr
    : std::conditional< 
	   traits::is_vector<Sub>, 
	   VectorExpr<Model>,
	   if_then_else< traits::is_matrix<Sub>, 
			 MatrixExpr<Model>,
			 BaseExpr<Model> > > {};

  
} // end namespace AMDiS
