/** \file MatrixVectorOperations.hpp */

#pragma once

#include "UnaryTerms.hpp"
#include "BinaryTerms.hpp"
#include "FunctorTerm.hpp"
#include "ConstantTerms.hpp"
#include "CoordsTerm.hpp"
#include "TermGenerator.hpp"

#include <traits/concepts.hpp>

namespace AMDiS 
{
  // ___________________________________________________________________________
  // generator functions for constants and references

  template <class T>
  constexpr RTConstant<T> constant(T&& value) { return {std::forward<T>(value)}; }
  
  template <int I>
  constexpr CTConstant<I> constant() { return {}; }
  
  template <class T>
  constexpr Reference<T> ref(T const& value) { return {value}; }
  
  inline CoordsTerm X() { return {}; } // can not be constexpr
  
  template <class F, class... Terms>
  constexpr requires::Term< FunctorTerm< F, ToTerm_t<Terms>... >, Terms... >
  func(F&& f, Terms&&... ts)
  {
    return {std::forward<F>(f), toTerm(std::forward<Terms>(ts))...}; 
  }
  
  // ---------------------------------------------------------------------------
  // Operations (elementwise)
  
  /// expression for V + W
  template <class T1, class T2>
  constexpr requires::Term< PlusTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2 > 
  operator+(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  
  /// expression for V - W
  template <class T1, class T2>
  constexpr requires::Term< MinusTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2 > 
  operator-(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  
  /// expression for -V
  template <class T>
  constexpr NegateTerm<T>
  operator-(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  
  /// expression for V .* W
  template <class T1, class T2>
  constexpr requires::Term< MultipliesTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2 >  
  operator*(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  
  /// expression for V ./ W
  template <class T1, class T2>
  constexpr requires::Term< DividesTerm<ToTerm<T1>, ToTerm<T2>>, T1, T2 > 
  operator/(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
#if 0

  // ---------------------------------------------------------------------------
  // scalar product
  
  /// scalar product V*V
  template <class T1, class T2>
  Value_t<typename DotExpr<T1, T2>::type>
  dot(VectorExpr<T1> const& term1, VectorExpr<T2> const& term2)
  {
    typedef typename DotExpr<T1, T2>::type binary_op;
    return binary_op(term1.sub(), term2.sub())();
  }
  
  /// expression for V * W (dot product)
  template <class T1, class T2>
  Value_t<typename DotExpr<T1, T2>::type>
  operator*(VectorExpr<T1> const& term1, VectorExpr<T2> const& term2)
  {
    return dot(term1, term2);
  }

  // T1 x T2
  template <class T1, class T2> // assume vector shape for T1 and T2
  struct CrossExpr {
    typedef VectorBinaryExpr<T1, T2, 
      functors::MyCross<Value_t<T1>, Value_t<T2>> > type;
  };
    
  /// expression for V x W (cross product / outer product / tensor product)
  template <class T1, class T2>
  typename CrossExpr<T1, T2>::type
  cross(VectorExpr<T1> const& term1, VectorExpr<T2> const& term2)
  {
    typedef typename CrossExpr<T1, T2>::type binary_op;
    return binary_op(term1.sub(), term2.sub());
  }
  
  // ---------------------------------------------------------------------------
  // reduction operations
  
  /// expression for one_norm(V)
  template <class E>
  Value_t<typename OneNormExpr<E>::type>
  one_norm(VectorExpr<E> const& expr)
  {
    typedef typename OneNormExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// expression for two_norm(V)
  template <class E>
  Value_t<typename TwoNormExpr<E>::type>
  two_norm(VectorExpr<E> const& expr)
  {
    typedef typename TwoNormExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// expression for two_norm(M)
  template <class E>
  Value_t<typename TwoNormExpr<E>::type>
  frobenius_norm(MatrixExpr<E> const& expr)
  {
    typedef typename TwoNormExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// expression for norm(V) := two_norm(V)
  template <class E>
  Value_t<typename TwoNormExpr<E>::type>
  norm(VectorExpr<E> const& expr)
  {
    return two_norm(expr);
  }
  
  /// expression for norm(M) := frobenius_norm(M)
  template <class E>
  Value_t<typename TwoNormExpr<E>::type>
  norm(MatrixExpr<E> const& expr)
  {
    return frobenius_norm(expr);
  }
  
  /// expression for unary_dot(V) = V*V
  template <class E>
  Value_t<typename UnaryDotExpr<E>::type>
  unary_dot(VectorExpr<E> const& expr)
  {
    typedef typename UnaryDotExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// expression for max(V)
  template <class E>
  Value_t<typename MaxExpr<E>::type>
  max(BaseTerm<E> const& expr)
  {
    typedef typename MaxExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// expression for abs_max(V)
  template <class E>
  Value_t<typename AbsMaxExpr<E>::type>
  abs_max(BaseTerm<E> const& expr)
  {
    typedef typename AbsMaxExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// expression for inf_norm(V) = abs_max(V)
  template <class E>
  Value_t<typename AbsMaxExpr<E>::type>
  inf_norm(VectorExpr<E> const& expr)
  {
    return abs_max(expr);
  }
  
  /// expression for inf_norm(M) = sqrt(rows*cols) * abs_max(M) = Gesamtnorm(M)
  template <class E>
  Value_t<typename AbsMaxExpr<E>::type>
  inf_norm(MatrixExpr<E> const& expr)
  {
    return abs_max(expr) * std::sqrt(size(expr));
  }
  
  /// expression for min(V)
  template <class E>
  Value_t<typename MinExpr<E>::type>
  min(BaseTerm<E> const& expr)
  {
    typedef typename MinExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// expression for abs_min(V)
  template <class E>
  Value_t<typename AbsMinExpr<E>::type>
  abs_min(BaseTerm<E> const& expr)
  {
    typedef typename AbsMinExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// expression for sum(V)
  template <class E>
  Value_t<typename SumExpr<E>::type>
  sum(BaseTerm<E> const& expr)
  {
    typedef typename SumExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// expression for sum(V) = v0 + v1 + v2 + ...
  template <class E>
  Value_t<typename SumExpr<E>::type>
  mean(BaseTerm<E> const& expr)
  {
    return sum(expr) / size(expr);
  }
  
  /// expression for prod(V) = v0 * v1 * v2 * ...
  template <class E>
  Value_t<typename ProdExpr<E>::type>
  prod(BaseTerm<E> const& expr)
  {
    typedef typename ProdExpr<E>::type op;
    return op(expr.sub())();
  }
  
  /// euklidean distance |V1 - V2|_2
  template <class T1, class T2>
  Value_t<typename UnaryDotExpr<T1>::type>
  distance(VectorExpr<T1> const& term1, VectorExpr<T2> const& term2) // NOTE: in AMDiS::absteukl
  {
    return std::sqrt(unary_dot(term1.sub() - term2.sub()));
  }
  
  
  // ---------------------------------------------------------------------------
  // matrix-vector multiplication
  
  /// expression for Mat * V
  template <class T1, class T2>
  MatVecExpr<T1, T2, false>
  operator*(MatrixExpr<T1> const& term1, VectorExpr<T2> const& term2)
  {
    typedef MatVecExpr<T1, T2, false> binary_op;
    return binary_op(term1.sub(), term2.sub());
  }
#endif
} // end namespace AMDiS
