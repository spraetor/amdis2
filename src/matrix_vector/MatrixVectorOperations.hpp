/** \file MatrixVectorOperations.hpp */

#pragma once

#include <matrix_vector/expr/all_expr.hpp>

#include <operations/functors.hpp>
#include <operations/reduction_functors.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>

#include <Math.h>

namespace AMDiS
{
  // ---------------------------------------------------------------------------
  // Operations with Vector and Matrix (elementwise)

  /// expression for V + W
  template <class E1, class E2>
  constexpr PlusExpr<E1, E2>
  operator+(BaseExpr<E1> const& expr1, BaseExpr<E2> const& expr2)
  {
    return {expr1.sub(), expr2.sub()};
  }


  /// expression for V - W
  template <class E1, class E2>
  constexpr MinusExpr<E1, E2>
  operator-(BaseExpr<E1> const& expr1, BaseExpr<E2> const& expr2)
  {
    return {expr1.sub(), expr2.sub()};
  }


  /// expression for -V
  template <class E>
  constexpr ElementwiseUnaryExpr<E, functors::negate<Value_t<E>>>
  operator-(BaseExpr<E> const& expr)
  {
    return {expr.sub()};
  }


  /// expression for V * scalar
  template <class Value, class E>
  constexpr RightScaleExpr<Value, E>
  operator*(BaseExpr<E> const& expr, Value scal)
  {
    return {scal, expr.sub()};
  }


  /// expression for scalar * V
  template <class Value, class E>
  constexpr LeftScaleExpr<Value, E>
  operator*(Value scal, BaseExpr<E> const& expr)
  {
    return {scal, expr.sub()};
  }


  /// expression for V / scalar
  template <class Value, class E>
  constexpr RightDivideExpr<Value, E>
  operator/(BaseExpr<E> const& expr, Value scal)
  {
    return {scal, expr.sub()};
  }

  // ---------------------------------------------------------------------------
  // scalar product

  /// scalar product V*V
  template <class E1, class E2>
  Value_t<DotExpr<E1, E2>>
  dot(VectorExpr<E1> const& expr1, VectorExpr<E2> const& expr2)
  {
    return DotExpr<E1, E2>(expr1.sub(), expr2.sub())();
  }

  /// expression for V * W (dot product)
  template <class E1, class E2>
  Value_t<DotExpr<E1, E2>>
  operator*(VectorExpr<E1> const& expr1, VectorExpr<E2> const& expr2)
  {
    return dot(expr1, expr2);
  }

  // E1 x E2
  template <class E1, class E2> // assume vector shape for E1 and E2
  using CrossExpr = VectorBinaryExpr<E1, E2,
			functors::MyCross<Value_t<E1>, Value_t<E2>>>;

  /// expression for V x W (cross product)
  template <class E1, class E2>
  CrossExpr<E1, E2>
  cross(VectorExpr<E1> const& expr1, VectorExpr<E2> const& expr2)
  {
    return {expr1.sub(), expr2.sub()};
  }

  // ---------------------------------------------------------------------------
  // reduction operations

  /// expression for one_norm(V)
  template <class E>
  Value_t<typename OneNormExpr<E>::type>
  one_norm(VectorExpr<E> const& expr)
  {
    using op = typename OneNormExpr<E>::type;
    return op(expr.sub())();
  }

  /// expression for two_norm(V)
  template <class E>
  Value_t<typename TwoNormExpr<E>::type>
  two_norm(VectorExpr<E> const& expr)
  {
    using op = typename TwoNormExpr<E>::type;
    return op(expr.sub())();
  }

  /// expression for two_norm(M)
  template <class E>
  Value_t<typename TwoNormExpr<E>::type>
  frobenius_norm(MatrixExpr<E> const& expr)
  {
    using op = typename TwoNormExpr<E>::type;
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
    using op = typename UnaryDotExpr<E>::type;
    return op(expr.sub())();
  }

  /// expression for max(V)
  template <class E>
  Value_t<typename MaxExpr<E>::type>
  max(BaseExpr<E> const& expr)
  {
    using op = typename MaxExpr<E>::type;
    return op(expr.sub())();
  }

  /// expression for abs_max(V)
  template <class E>
  Value_t<typename AbsMaxExpr<E>::type>
  abs_max(BaseExpr<E> const& expr)
  {
    using op = typename AbsMaxExpr<E>::type;
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
    return abs_max(expr) * std::sqrt(static_cast<double>(size(expr)));
  }

  /// expression for min(V)
  template <class E>
  Value_t<typename MinExpr<E>::type>
  min(BaseExpr<E> const& expr)
  {
    using op = typename MinExpr<E>::type;
    return op(expr.sub())();
  }

  /// expression for abs_min(V)
  template <class E>
  Value_t<typename AbsMinExpr<E>::type>
  abs_min(BaseExpr<E> const& expr)
  {
    using op = typename AbsMinExpr<E>::type;
    return op(expr.sub())();
  }

  /// expression for sum(V)
  template <class E>
  Value_t<typename SumExpr<E>::type>
  sum(BaseExpr<E> const& expr)
  {
    using op = typename SumExpr<E>::type;
    return op(expr.sub())();
  }

  /// expression for sum(V) = v0 + v1 + v2 + ...
  template <class E>
  Value_t<typename SumExpr<E>::type>
  mean(BaseExpr<E> const& expr)
  {
    return sum(expr) / size(expr);
  }

  /// expression for prod(V) = v0 * v1 * v2 * ...
  template <class E>
  Value_t<typename ProdExpr<E>::type>
  prod(BaseExpr<E> const& expr)
  {
    using op = typename ProdExpr<E>::type;
    return op(expr.sub())();
  }

  /// euklidean distance |V1 - V2|_2
  template <class E1, class E2>
  Value_t<typename UnaryDotExpr<E1>::type>
  distance(VectorExpr<E1> const& expr1, VectorExpr<E2> const& expr2)
  {
    return std::sqrt(unary_dot(expr1.sub() - expr2.sub()));
  }


  // ---------------------------------------------------------------------------
  // matrix-vector multiplication

  /// expression for Mat * V
  template <class E1, class E2>
  MatVecExpr<E1, E2, false>
  operator*(MatrixExpr<E1> const& expr1, VectorExpr<E2> const& expr2)
  {
    return {expr1.sub(), expr2.sub()};
  }


  /// comparison of expressions
  template <class E1, class E2>
  inline bool operator==(VectorExpr<E1> const& expr1, VectorExpr<E2> const& expr2)
  {
    E1 const& v1 = expr1.sub();
    E2 const& v2 = expr2.sub();

    for (Size_t<E1> i = 0; i < size(v1); ++i)
      if (math::abs(v1(i) - v2(i)) > DBL_TOL)
        return false;

    return false;
  }


  /// comparison of expressions
  template <class E1, class E2>
  inline bool operator!=(VectorExpr<E1> const& expr1, VectorExpr<E2> const& expr2)
  {
    return !(expr1 == expr2);
  }


  /// test for less-then (elementwise) up to DBL_TOL
  template <class E1, class E2>
  inline bool operator<(VectorExpr<E1> const& expr1, VectorExpr<E2> const& expr2)
  {
    E1 const& v1 = expr1.sub();
    E2 const& v2 = expr2.sub();

    for (Size_t<E1> i = 0; i < size(v1); ++i)
    {
      if (math::abs(v1(i) - v2(i)) < DBL_TOL)
        continue;
      return v1(i) < v2(i);
    }
    return false;
  }
} // end namespace AMDiS
