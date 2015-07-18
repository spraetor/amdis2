/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file vec_functors.hpp */

#ifndef AMDIS_VEC_FUNCTORS_HPP
#define AMDIS_VEC_FUNCTORS_HPP

#include "functor_expr.hpp"
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
    template<typename T> struct is_always_true : boost::mpl::true_ {};
  }
  
  namespace result_of
  {
    // helper class for UnaryExpr
    template<template<class> class Functor,
	     typename Term,
	     template<class> class Condition1 = traits::is_always_true,
	     template<class> class Condition2 = traits::is_always_true,
	     template<class> class Condition3 = traits::is_always_true>
    struct UnaryExpr : boost::enable_if
      <
	typename boost::mpl::and_
	<
	  typename traits::is_expr<Term>::type,
	  typename boost::mpl::or_
	  <
	    typename Condition1<typename Term::value_type>::type,
	    typename Condition2<typename Term::value_type>::type,
	    typename Condition3<typename Term::value_type>::type
	  >::type
	>,
	expressions::Function1<Functor<typename Term::value_type>, Term>
      > {};
      
    // helper class for UnaryExpr
    template<class Functor,
	     typename Term,
	     template<class> class Condition1 = traits::is_always_true,
	     template<class> class Condition2 = traits::is_always_true,
	     template<class> class Condition3 = traits::is_always_true>
    struct UnaryExprFull : boost::enable_if
      <
	typename boost::mpl::and_
	<
	  typename traits::is_expr<Term>::type,
	  typename boost::mpl::or_
	  <
	    typename Condition1<typename Term::value_type>::type,
	    typename Condition2<typename Term::value_type>::type,
	    typename Condition3<typename Term::value_type>::type
	  >::type
	>,
	expressions::Function1<Functor, Term>
      > {};
      
      
    // helper class for BinaryExpr
    template<template<class,class> class Functor,
	     typename Term1, typename Term2, 
	     template<class> class Condition1 = traits::is_always_true,
	     template<class> class Condition2 = traits::is_always_true,
	     template<class> class Condition3 = traits::is_always_true,
	     typename Expr1 = typename traits::to_expr<Term1>::type,
	     typename Expr2 = typename traits::to_expr<Term2>::type>
    struct BinaryExpr : boost::enable_if
      < 
	typename boost::mpl::and_
	<
	  typename traits::is_valid_arg2<Term1, Term2>::type,
	  typename boost::mpl::or_
	  <
	    typename Condition1<typename Expr1::value_type>::type,
	    typename Condition2<typename Expr1::value_type>::type,
	    typename Condition3<typename Expr1::value_type>::type
	  >::type,
	  typename boost::mpl::or_
	  <
	    typename Condition1<typename Expr2::value_type>::type,
	    typename Condition2<typename Expr2::value_type>::type,
	    typename Condition3<typename Expr2::value_type>::type
	  >::type
	>::type,
	expressions::Function2
	<
	  Functor<typename Expr1::value_type, typename Expr2::value_type>,
	  Expr1, Expr2
	> 
      > {};
      
  } // end namespace result_of
  
    
  /// expression with one vector
  template<template<class> class Functor, typename Term>
  inline typename result_of::UnaryExpr<Functor, Term>::type
  unary_expr(const Term& t) 
  { 
    typedef Functor< typename Term::value_type > F;
    return function_(F(), t); 
  }
  
  template<class F, typename Term>
  inline typename result_of::UnaryExprFull<F, Term>::type
  unary_expr_full(const Term& t) 
  { 
    return function_(F(), t); 
  }
  
    
  /// expression with two vectors v1, v2
  template<template<class,class> class Functor, typename Term1, typename Term2>
  inline typename result_of::BinaryExpr<Functor, Term1, Term2>::type
  binary_expr(const Term1& t1, const Term2& t2) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    typedef Functor< typename Expr1::type::value_type, typename Expr2::type::value_type > F;
    return function_(F(), Expr1::get(t1), Expr2::get(t2)); 
  }
    

// two_norm of vectors
// _____________________________________________________________________________

  namespace result_of 
  {
    template<typename T>
    struct two_norm<T, tag::expression> : result_of::UnaryExpr<functors::TwoNorm, T, traits::is_vector> {};
  }
  
  /// the 2-norm of a vector v: result = sqrt(v^H * v)
  template<typename Term>
  inline typename result_of::two_norm<Term, tag::expression>::type
  two_norm_dispatch(const Term& t, tag::expression)
  {
    return unary_expr< functors::TwoNorm >(t); 
  }

  
// one_norm of vectors and matrices
// _____________________________________________________________________________

  namespace result_of 
  {
    template<typename T>
    struct one_norm<T, tag::expression> : result_of::UnaryExpr<functors::OneNorm, T, traits::is_vector, traits::is_matrix> {};
  }
  
  /// the 1-norm of a vector v: result = max_i(abs(v_i))
  template<typename Term>
  inline typename result_of::one_norm<Term, tag::expression>::type
  one_norm_dispatch(const Term& t, tag::expression)
  {
    return unary_expr< functors::OneNorm >(t); 
  }

  
// p_norm<P> of vectors
// _____________________________________________________________________________

  namespace result_of 
  {
    template<int P, typename T>
    struct p_norm<P, T, tag::expression> : result_of::UnaryExprFull<functors::PNorm<P, typename T::value_type>, T, traits::is_vector> {};
  }
  
  /// the 2-norm of a vector v: result = sqrt(v^H * v)
  template<int P, typename Term>
  inline typename result_of::p_norm<P, Term, tag::expression>::type
  p_norm_dispatch(const Term& t, tag::expression)
  {
    return unary_expr_full< functors::PNorm<P, typename Term::value_type> >(t); 
  }

  
// unary dot (v' * v) of a vector v
// _____________________________________________________________________________

  namespace result_of 
  {
    template<typename T>
    struct unary_dot<T, tag::expression> : result_of::UnaryExpr<functors::UnaryDot, T, traits::is_vector> {};
  }
  
  /// the 2-norm squared of a vector v: result = v^H * v
  template<typename Term>
  inline typename result_of::unary_dot<Term, tag::expression>::type
  unary_dot_dispatch(const Term& t, tag::expression)
  {
    return unary_expr< functors::UnaryDot >(t); 
  }
  

// inner product of two vectors
// _____________________________________________________________________________
  
  namespace result_of 
  {
    template<typename T1, typename T2>
    struct dot<T1, T2, tag::expression, tag::expression> : result_of::BinaryExpr<functors::Inner, T1, T2, traits::is_vector> {};
    template<typename T1, typename T2>
    struct dot<T1, T2, tag::expression, tag::vector> : result_of::BinaryExpr<functors::Inner, T1, T2, traits::is_vector> {};
    template<typename T1, typename T2>
    struct dot<T1, T2, tag::vector, tag::expression> : result_of::BinaryExpr<functors::Inner, T1, T2, traits::is_vector> {};
  }
  
  /// inner/dot product of two vectors v1, v2: result = v1^T * v2
  template<typename Term1, typename Term2, typename Tag2>
  inline typename boost::enable_if< typename traits::is_expr<Term1>::type,
    typename result_of::dot<Term1, Term2, tag::expression, Tag2>::type
    >::type
  dot_dispatch(const Term1& t1, const Term2& t2, tag::expression tag1_, Tag2 tag2_)
  {
    return binary_expr<functors::Inner>(t1, t2);
  }
  
  template<typename Term1, typename Term2, typename Tag1>
  inline typename boost::enable_if_c< (traits::is_expr<Term2>::type::value && !(traits::is_expr<Term1>::type::value)),
    typename result_of::dot<Term1, Term2, Tag1, tag::expression>::type
    >::type
  dot_dispatch(const Term1& t1, const Term2& t2, Tag1 tag1_, tag::expression tag2_)
  {
    return binary_expr<functors::Inner>(t1, t2);
  }

  
// cross-product of two vectors
// _____________________________________________________________________________
  
  namespace result_of 
  {
    template<typename T1, typename T2>
    struct cross<T1, T2, tag::expression, tag::expression> : result_of::BinaryExpr<functors::Cross, T1, T2, traits::is_vector> {};
    template<typename T1, typename T2>
    struct cross<T1, T2, tag::expression, tag::vector> : result_of::BinaryExpr<functors::Cross, T1, T2, traits::is_vector> {};
    template<typename T1, typename T2>
    struct cross<T1, T2, tag::vector, tag::expression> : result_of::BinaryExpr<functors::Cross, T1, T2, traits::is_vector> {};
  }

  /// cross product of two vectors v1, v2: result = v1 x v2     
  template<typename Term1, typename Term2, typename Tag2>
  inline typename boost::enable_if< typename traits::is_expr<Term1>::type,
    typename result_of::cross<Term1, Term2, tag::expression, Tag2>::type
    >::type
  cross_dispatch(const Term1& t1, const Term2& t2, tag::expression tag1_, Tag2 tag2_)
  {
    return binary_expr<functors::Cross>(t1, t2);
  }
  
  template<typename Term1, typename Term2, typename Tag1>
  inline typename boost::enable_if_c< (traits::is_expr<Term2>::type::value && !(traits::is_expr<Term1>::type::value)),
    typename result_of::cross<Term1, Term2, Tag1, tag::expression>::type
    >::type
  cross_dispatch(const Term1& t1, const Term2& t2, Tag1 tag1_, tag::expression tag2_)
  {
    return binary_expr<functors::Cross>(t1, t2);
  }

  
// generator for a diagonal matrix from a vector
// _____________________________________________________________________________

  namespace expressions 
  {
    template<typename VectorType>
    struct DiagonalMat : public FunctorBase {};
    
    template<typename T>
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
    
    template<typename T>
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
    
    template<typename T>
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
    
    template<typename T, typename Param>
    struct DiagonalMat<mtl::dense_vector<T, Param> > : public FunctorBase
    {
      typedef mtl::matrix::compressed2D<T, mtl::matrix::parameters<> > result_type;
      
      int getDegree(int d0) const { return d0; }
      template<typename Vector>
      result_type operator()(const Vector &v) const 
      { 
	return diagonal(v); 
      }
    };
  }


  /// create diagonal matrix from vector
  template<typename Term>
  inline typename result_of::UnaryExpr<expressions::DiagonalMat, Term, traits::is_vector>::type
  diagonal(const Term& t)
  { 
    return unary_expr<expressions::DiagonalMat>(t); 
  }

  
// outer product / dyadic product of two vectors to create a matrix
// _____________________________________________________________________________
  
  /// outer product of two vectors v1, v2: result = v1 * v2^T
  template<typename Term1, typename Term2>
  inline typename result_of::BinaryExpr<functors::Outer, Term1, Term2, traits::is_vector>::type
  outer(const Term1& t1, const Term2& t2) 
  {
    return binary_expr<functors::Outer>(t1, t2);
  }

  
// extract a component of a vector/matrix
// _____________________________________________________________________________

  namespace expressions 
  {
    template<typename T>
    struct VecComponent : public FunctorBase
    {
      typedef typename traits::category<T>::value_type result_type;
      VecComponent(int I_) : I(I_) {}
      
      int getDegree(int d0) const { return d0; }
      result_type operator()(const T &v) const { return v[I]; }
    private:
      int I;
    };
    
    template<typename T>
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

  template<typename Term>
  typename result_of::UnaryExpr<expressions::VecComponent, Term, traits::is_vector>::type
  at(const Term& t, int I) 
  {
    return function_(expressions::VecComponent<typename Term::value_type>(I), t); 
  }

  template<typename Term>
  typename result_of::UnaryExpr<expressions::MatComponent, Term, traits::is_matrix>::type
  at(const Term& t, int I, int J) 
  {
    return function_(expressions::MatComponent<typename Term::value_type>(I, J), t); 
  }


  
// transpose a matrix
// _____________________________________________________________________________

  namespace expressions 
  {
    template<typename Mat>
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
  }

  template<typename Term>
  typename result_of::UnaryExpr<expressions::MatTranspose, Term, traits::is_matrix>::type
  trans(const Term& t) 
  {
    return function_(expressions::MatTranspose<typename Term::value_type>(), t); 
  }
}

#endif // AMDIS_VEC_FUNCTORS_HPP
