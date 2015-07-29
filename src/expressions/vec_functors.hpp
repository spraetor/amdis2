/** \file vec_functors.hpp */

#pragma once

#include "functorN_expr.hpp"
#include "operations/norm.hpp"
#include "operations/product.hpp"
    
/**
 *  This file provides expressions for vectors and matrices:
 *  (where v_expr is an expression of vector type, and m_expr an
 *   expression of matrix type)
 *
 *    two_norm(v_expr) ... the 2-norm of a vector v: result = sqrt(v^H * v)
 *    one_norm(v_expr) ... the 1-norm of a vector v: result = sum_i(abs(v_i))
 *    one_norm(m_expr) ... the 1-norm of a matrix m: result = max_j(sum_i(abs(m_ij))
 *    p_norm<P>(v_expr) .. the P-norm of a vector v: result = [sum_i(abs(v_i)^P)]^(1/P)
 *
 **/    
    
namespace AMDiS 
{
  namespace traits
  {
    template <class T> using is_always_true = true_;
  }
  
  namespace result_of
  {
    // helper class for UnaryExpr
    template<template<class> class Functor,
	     class Term,
	     template<class> class Condition1 = traits::is_always_true,
	     template<class> class Condition2 = traits::is_always_true,
	     template<class> class Condition3 = traits::is_always_true>
    using UnaryExpr = enable_if
      <
      	and_<
      	  traits::is_expr<Term>,
      	  or_< Condition1<Value_t<Term> >,
      	       Condition2<Value_t<Term> >,
      	       Condition3<Value_t<Term> > >
      	>,
      	expressions::Function1<Functor<Value_t<Term> >, Term>
      >;
      
    // helper class for UnaryExpr
    template<class Functor,
	     class Term,
	     template<class> class Condition1 = traits::is_always_true,
	     template<class> class Condition2 = traits::is_always_true,
	     template<class> class Condition3 = traits::is_always_true>
    using UnaryExprFull = enable_if
      <
      	and_<
      	  traits::is_expr<Term>,
      	  or_< Condition1<Value_t<Term> >,
      	       Condition2<Value_t<Term> >,
      	       Condition3<Value_t<Term> > >
      	>,
      	expressions::Function1<Functor, Term>
      >;
      
      
    // helper class for BinaryExpr
    template<template<class,class> class Functor,
	     class Term1, class Term2, 
	     template<class> class Condition1 = traits::is_always_true,
	     template<class> class Condition2 = traits::is_always_true,
	     template<class> class Condition3 = traits::is_always_true,
	     class Expr1 = typename traits::to_expr<Term1>::type,
	     class Expr2 = typename traits::to_expr<Term2>::type>
    using BinaryExpr = enable_if
      < 
      	and_<
      	  traits::is_valid_arg2<Term1, Term2>,
      	  or_< Condition1<Value_t<Expr1> >,
      	       Condition2<Value_t<Expr1> >,
      	       Condition3<Value_t<Expr1> > >,
               
      	  or_< Condition1<Value_t<Expr2> >,
      	       Condition2<Value_t<Expr2> >,
      	       Condition3<Value_t<Expr2> > >
      	>,
      	expressions::Function2
      	<
      	  Functor<Value_t<Expr1>, Value_t<Expr2> >,
      	  Expr1, Expr2
      	> 
      >;
      
  } // end namespace result_of
  
    
  /// expression with one vector
  template <template<class> class Functor, class Term>
  inline typename result_of::UnaryExpr<Functor, Term>::type
  unary_expr(Term&& t) 
  { 
    using F = Functor< Value_t<Term> >;
    return function_(F(), std::forward<Term>(t)); 
  }
  
  template <class F, class Term>
  inline typename result_of::UnaryExprFull<F, Term>::type
  unary_expr_full(Term&& t) 
  { 
    return function_(F(), std::forward<Term>(t)); 
  }
  
    
  /// expression with two vectors v1, v2
  template <template<class,class> class Functor, class Term1, class Term2>
  inline typename result_of::BinaryExpr<Functor, Term1, Term2>::type
  binary_expr(Term1&& t1, Term2&& t2) 
  { 
    using Expr1 = traits::to_expr<Term1>;
    using Expr2 = traits::to_expr<Term2>;
    using F = Functor< Value_t<typename Expr1::type>, Value_t<typename Expr2::type> >;
    return function_(F(), Expr1::get(t1), Expr2::get(t2)); 
  }
    

// two_norm of vectors
// _____________________________________________________________________________

  namespace result_of 
  {
    template <class T>
    struct two_norm<T, tag::expression> 
      : result_of::UnaryExpr<functors::TwoNorm, T, traits::is_vector> {};
  }
  
  /// the 2-norm of a vector v: result = sqrt(v^H * v)
  template <class Term>
  inline typename result_of::two_norm<Term, tag::expression>::type
  two_norm_dispatch(Term&& t, tag::expression)
  {
    return unary_expr< functors::TwoNorm >(std::forward<Term>(t)); 
  }

  
// one_norm of vectors and matrices
// _____________________________________________________________________________

  namespace result_of 
  {
    template <class T>
    struct one_norm<T, tag::expression> 
      : result_of::UnaryExpr<functors::OneNorm, T, traits::is_vector, traits::is_matrix> {};
  }
  
  /// the 1-norm of a vector v: result = max_i(abs(v_i))
  template <class Term>
  inline typename result_of::one_norm<Term, tag::expression>::type
  one_norm_dispatch(Term&& t, tag::expression)
  {
    return unary_expr< functors::OneNorm >(std::forward<Term>(t)); 
  }

  
// p_norm<P> of vectors
// _____________________________________________________________________________

  namespace result_of 
  {
    template <int P, class T>
    struct p_norm<P, T, tag::expression> 
      : result_of::UnaryExprFull<functors::PNorm<P, Value_t<T> >, T, traits::is_vector> {};
  }
  
  /// the 2-norm of a vector v: result = sqrt(v^H * v)
  template <int P, class Term>
  inline typename result_of::p_norm<P, Term, tag::expression>::type
  p_norm_dispatch(Term&& t, tag::expression)
  {
    return unary_expr_full< functors::PNorm<P, Value_t<Term> > >(std::forward<Term>(t)); 
  }

  
// unary dot (v' * v) of a vector v
// _____________________________________________________________________________

  namespace result_of 
  {
    template <class T>
    struct unary_dot<T, tag::expression> 
      : result_of::UnaryExpr<functors::UnaryDot, T, traits::is_vector> {};
  }
  
  /// the 2-norm squared of a vector v: result = v^H * v
  template <class Term>
  inline typename result_of::unary_dot<Term, tag::expression>::type
  unary_dot_dispatch(Term&& t, tag::expression)
  {
    return unary_expr< functors::UnaryDot >(std::forward<Term>(t)); 
  }
  

// inner product of two vectors
// _____________________________________________________________________________
  
  namespace result_of 
  {
    template <class T1, class T2>
    struct dot<T1, T2, tag::expression, tag::expression> 
      : result_of::BinaryExpr<functors::Inner, T1, T2, traits::is_vector> {};
    template <class T1, class T2>
    struct dot<T1, T2, tag::expression, tag::vector> 
      : result_of::BinaryExpr<functors::Inner, T1, T2, traits::is_vector> {};
    template <class T1, class T2>
    struct dot<T1, T2, tag::vector, tag::expression> 
      : result_of::BinaryExpr<functors::Inner, T1, T2, traits::is_vector> {};
  }
  
  /// inner/dot product of two vectors v1, v2: result = v1^T * v2
  template <class Term1, class Term2, class Tag2>
  inline typename enable_if< traits::is_expr<Term1>,
    typename result_of::dot<Term1, Term2, tag::expression, Tag2>::type
    >::type
  dot_dispatch(Term1&& t1, Term2&& t2, tag::expression tag1_, Tag2 tag2_)
  {
    return binary_expr<functors::Inner>(std::forward<Term1>(t1), std::forward<Term2>(t2));
  }
  
  template <class Term1, class Term2, class Tag1>
  inline typename enable_if_c< (traits::is_expr<Term2>::value && !(traits::is_expr<Term1>::value)),
    typename result_of::dot<Term1, Term2, Tag1, tag::expression>::type
    >::type
  dot_dispatch(Term1&& t1, Term2&& t2, Tag1 tag1_, tag::expression tag2_)
  {
    return binary_expr<functors::Inner>(std::forward<Term1>(t1), std::forward<Term2>(t2));
  }

  
// cross-product of two vectors
// _____________________________________________________________________________
  
  namespace result_of 
  {
    template<class T1, class T2>
    struct cross<T1, T2, tag::expression, tag::expression> 
      : result_of::BinaryExpr<functors::Cross, T1, T2, traits::is_vector> {};
    template<class T1, class T2>
    struct cross<T1, T2, tag::expression, tag::vector> 
      : result_of::BinaryExpr<functors::Cross, T1, T2, traits::is_vector> {};
    template<class T1, class T2>
    struct cross<T1, T2, tag::vector, tag::expression> 
      : result_of::BinaryExpr<functors::Cross, T1, T2, traits::is_vector> {};
  }

  /// cross product of two vectors v1, v2: result = v1 x v2     
  template <class Term1, class Term2, class Tag2>
  inline typename enable_if< traits::is_expr<Term1>,
    typename result_of::cross<Term1, Term2, tag::expression, Tag2>::type
    >::type
  cross_dispatch(Term1&& t1, Term2&& t2, tag::expression tag1_, Tag2 tag2_)
  {
    return binary_expr<functors::Cross>(std::forward<Term1>(t1), std::forward<Term2>(t2));
  }
  
  template <class Term1, class Term2, class Tag1>
  inline typename enable_if_c< (traits::is_expr<Term2>::value && !(traits::is_expr<Term1>::value)),
    typename result_of::cross<Term1, Term2, Tag1, tag::expression>::type
    >::type
  cross_dispatch(Term1&& t1, Term2&& t2, Tag1 tag1_, tag::expression tag2_)
  {
    return binary_expr<functors::Cross>(std::forward<Term1>(t1), std::forward<Term2>(t2));
  }

  
// generator for a diagonal matrix from a vector
// _____________________________________________________________________________

  namespace expressions 
  {
    template <class VectorType>
    struct DiagonalMat : public FunctorBase {};
    
    template <class T>
    struct DiagonalMat<WorldVector<T> > : public FunctorBase
    {
      typedef WorldMatrix<T> result_type;
      
      int getDegree(int d0) const { return d0; }
      result_type operator()(const WorldVector<T> &v) const 
      { 
      	T zero; nullify(zero);
      	result_type matrix(DEFAULT_VALUE, zero);
      	for (int i = 0; i < v.getSize(); ++i)
      	  matrix[i][i] = v[i];
      	return matrix; 
      }
    };
    
    template <class T>
    struct DiagonalMat<DimVec<T> > : public FunctorBase
    {
      typedef DimMat<T> result_type;
      
      int getDegree(int d0) const { return d0; }
      result_type operator()(const DimVec<T> &v) const 
      { 
      	T zero; nullify(zero);
      	result_type matrix(v.getDim(), DEFAULT_VALUE, zero);
      	for (int i = 0; i < v.getSize(); ++i)
      	  matrix[i][i] = v[i];
      	return matrix; 
      }
    };
    
    template <class T>
    struct DiagonalMat<Vector<T> > : public FunctorBase
    {
      typedef Matrix<T> result_type;
      
      int getDegree(int d0) const { return d0; }
      result_type operator()(const DimVec<T> &v) const 
      { 
      	T zero; nullify(zero);
      	result_type matrix(v.getSize(), v.getSize());
      	matrix.set(zero);
      	for (int i = 0; i < v.getSize(); ++i)
      	  matrix[i][i] = v[i];
      	return matrix; 
      }
    };
    
    template <class T, class Param>
    struct DiagonalMat<mtl::dense_vector<T, Param> > : public FunctorBase
    {
      typedef mtl::matrix::compressed2D<T, mtl::matrix::parameters<> > result_type;
      
      int getDegree(int d0) const { return d0; }
      template <class Vector>
      result_type operator()(const Vector &v) const 
      { 
        return diagonal(v); 
      }
    };
  }


  /// create diagonal matrix from vector
  template <class Term>
  inline typename result_of::UnaryExpr<expressions::DiagonalMat, Term, traits::is_vector>::type
  diagonal(Term&& t)
  { 
    return unary_expr<expressions::DiagonalMat>(std::forward<Term>(t)); 
  }

  
// outer product / dyadic product of two vectors to create a matrix
// _____________________________________________________________________________
  
  /// outer product of two vectors v1, v2: result = v1 * v2^T
  template <class Term1, class Term2>
  inline typename result_of::BinaryExpr<functors::Outer, Term1, Term2, traits::is_vector>::type
  outer(Term1&& t1, Term2&& t2) 
  {
    return binary_expr<functors::Outer>(std::forward<Term1>(t1), std::forward<Term2>(t2));
  }

  
// extract a component of a vector/matrix
// _____________________________________________________________________________

  namespace expressions 
  {
    template <class T>
    struct VecComponent : public FunctorBase
    {
      typedef typename traits::category<T>::value_type result_type;
      VecComponent(int I_) : I(I_) {}
      
      int getDegree(int d0) const { return d0; }
      result_type operator()(const T &v) const { return v[I]; }
    private:
      int I;
    };
    
    template <class T>
    struct MatComponent : public FunctorBase
    {
      typedef typename traits::category<T>::value_type result_type;
      MatComponent(int I_, int J_) : I(I_), J(J_) {}
      
      int getDegree(int d0) const { return d0; }
      result_type operator()(const T &m) const { return m[I][J]; }
    private:
      int I;
      int J;
    };
  }

  template <class Term>
  typename result_of::UnaryExpr<expressions::VecComponent, Term, traits::is_vector>::type
  at(Term&& t, int I) 
  {
    return function_(expressions::VecComponent<Value_t<Term>>(I), std::forward<Term>(t)); 
  }

  template <class Term>
  typename result_of::UnaryExpr<expressions::MatComponent, Term, traits::is_matrix>::type
  at(Term&& t, int I, int J) 
  {
    return function_(expressions::MatComponent<Value_t<Term>>(I, J), std::forward<Term>(t)); 
  }


  
// transpose a matrix
// _____________________________________________________________________________

  namespace expressions 
  {
    template <class Mat>
    struct MatTranspose : public FunctorBase
    {
      typedef Mat result_type;
      int getDegree(int d0) const { return d0; }
      result_type operator()(const Mat &m) const 
      { 
        Mat result;
        for (size_t r = 0; r < num_rows(m); r++)
          for (size_t c = 0; c < num_cols(m); c++)
            result[c][r] = m[r][c];
        return result; 
      }
    };
    
  } // end namespace expressions

  template<typename Term>
  typename result_of::UnaryExpr<expressions::MatTranspose, Term, traits::is_matrix>::type
  trans(Term&& t) 
  {
    return function_(expressions::MatTranspose<Value_t<Term>>(), std::forward<Term>(t)); 
  }
  
} // end namespace AMDiS
