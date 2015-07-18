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



/** \file MatrixVectorOperations.h */

#ifndef AMDIS_MATVEC_OPERATIONS_H
#define AMDIS_MATVEC_OPERATIONS_H

#include "Traits.h"

namespace AMDiS {

  // ---------------------------------------------------------------------------
  // Operations with Vector and Matrix

  /// Matrix vector multiplication.
  template<typename T>
  inline const Vector<T>& mv(const Matrix<T>& m, const Vector<T>& v, Vector<T>& result)
  {
    TEST_EXIT_DBG(m.getNumCols() == v.getSize())("m and v not compatible\n");
    TEST_EXIT_DBG(v.getSize() == result.getSize())("wrong result size\n");

    T *resultIt;
    T const* mIt;
    T const* vIt;

    for (resultIt = result.begin(), mIt = m.begin(); 
	 resultIt != result.end(); 
	 ++resultIt) {
      *resultIt = 0;
      for (vIt = v.begin(); vIt != v.end(); ++vIt, ++mIt)
	*resultIt += *mIt * *vIt;
    }

    return result;
  }

  /// Vector addition.
  template<typename T> 
  inline const Vector<T>& add(const Vector<T>& v1, const Vector<T>& v2, Vector<T>& result)
  {
    TEST_EXIT_DBG(v1.getSize() == v2.getSize())("invalid size in test v1 == v2\n");
    TEST_EXIT_DBG(v2.getSize() == result.getSize())("invalid size in test v2 == result\n");
    T const* v1It;
    T const* v2It;
    T* resultIt;
    for (v1It = v1.begin(), v2It = v2.begin(), resultIt = result.begin();
	 v1It != v1.end();
	 ++v1It, ++v2It, ++resultIt)
      *resultIt = *v1It + *v2It;

    return result;
  }

  /// scalar * vector
  template<typename T, typename S>
  inline const Vector<T>& mult(const S& scal,
			       const Vector<T>& v, 
			       Vector<T>& result)
  {
    TEST_EXIT_DBG(v.getSize() == result.getSize())("invalid size\n");

    T const* vIt;
    T* resultIt;
    for (vIt = v.begin(), resultIt = result.begin();
	 vIt != v.end();
	 ++vIt, ++resultIt) 
      *resultIt = scal * *vIt;

    return result;
  }

  /// vector + scalar
  template<typename T>
  inline const Vector<T>& add(const Vector<T>& v, const T& scal, Vector<T>& result)
  {
    TEST_EXIT_DBG(v.getSize() == result.getSize())("invalid size\n");
    T const* vIt;
    T* resultIt;
    for (vIt = v.begin(), resultIt = result.begin();
	 vIt != v.end();
	 ++vIt, ++resultIt)
      *resultIt = *vIt + scal;

    return result;
  }

  /// y = a * x + y.
  template<typename T, typename S>
  inline const Vector<T>& axpy(const S& a,
			       const Vector<T> &x,
			       Vector<T> &y)
  {
    TEST_EXIT_DBG(x.getSize() == y.getSize())("invalid size\n");
    T const* xIt;
    T* yIt;
    for (xIt = x.begin(), yIt = y.begin();
	 xIt != x.end();
	 ++xIt, ++yIt) 
      *yIt += a * *xIt;

    return y;
  }
  
  // times / divides
  // ---------------

  /// Matrix vector multiplication: vector := matrix * vector
  template<typename T>
  inline const Vector<T>& operator*=(const Vector<T>& v, const Matrix<T>& m) 
  {
    return mv(m, v, v);
  }

  /// Matrix vector multiplication: vector := matrix * vector
  template<typename T>
  inline Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v) 
  {
    Vector<T> result(m.getNumCols());
    return mv(m, v, result);
  }

  /// Scalar product: scalar := vector * vector
  template<typename T, typename S> 
  inline typename traits::mult_type<T,S>::type 
  operator*(const Vector<T>& v1, const Vector<S>& v2) 
  {
    typename traits::mult_type<T,S>::type result;
    nullify(result);

    T const* v1It;
    S const* v2It;
    for (v1It = v1.begin(), v2It = v2.begin();
	 v1It != v1.end();
	 ++v1It, ++v2It)
      result += *v1It * *v2It;

    return result;
  }

  /// vector *= scalar (elementwise)
  template <typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, Vector<T> >::type &
  operator*=(Vector<T>& v, S scal)
  {
    T *vIt;
    for (vIt = v.begin(); vIt != v.end(); ++vIt) 
      *vIt *= scal;
    return v;
  }

  /// vector := vector * scalar (elementwise)
  template <typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, Vector<T> >::type
  operator*(Vector<T> result, S scal) 
  {
    result *= scal;
    return result;
  }

  /// vector /= scalar (elementwise)
  template <typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, Vector<T> >::type &
  operator/=(Vector<T>& v, S scal)
  {
    T *vIt;
    for (vIt = v.begin(); vIt != v.end(); ++vIt) 
      *vIt /= scal;
    return v;
  }

  /// vector := vector / scalar (elementwise)
  template <typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, Vector<T> >::type
  operator/(Vector<T> result, S scal) 
  {
    result /= scal;
    return result;
  }
  
  // plus / minus
  // ------------

  /// vector += scalar
  template <typename T, typename S>
  typename boost::enable_if<typename traits::is_scalar<S>::type,
    Vector<T> >::type&
  operator+=(Vector<T>& x, S scal) 
  {
    T *xIt;
    for (xIt = x.begin(); xIt != x.end(); ++xIt) 
      *xIt += scal;
    return x;
  }

  /// vector := vector + scalar
  template <typename T, typename S>
  typename boost::enable_if<typename traits::is_scalar<S>::type,
    Vector<T> >::type
  operator+(Vector<T> result, T scal) 
  {
    result += scal;
    return result;
  }

  /// vector += vector
  template <typename T, typename S>
  Vector<T>& operator+=(Vector<T>& x, const Vector<S>& y) 
  {
    T *xIt;
    S const* yIt;
    for (xIt = x.begin(), yIt = y.begin(); xIt != x.end(); ++xIt, ++yIt) 
      *xIt += *yIt;
    return x;
  }

  /// vector := vector + vector
  template <typename T, typename S>
  Vector<T> operator+(Vector<T> result, const Vector<S>& v2) 
  {
    result += v2;
    return result;
  }

  /// vector -= vector
  template <typename T, typename S>
  Vector<T>& operator-=(Vector<T>& x, const Vector<S>& y)
  {
    T* xIt;
    S const* yIt;
    for (xIt = x.begin(), yIt = y.begin(); xIt != x.end(); ++xIt, ++yIt) 
      *xIt -= *yIt;
    return x;
  }

  /// vector := vector - vector
  template <typename T, typename S>
  Vector<T> operator-(Vector<T> result, const Vector<S>& v2)
  {
    result -= v2;
    return result;
  }
  
  // special operators
  // -----------------

  /// 2-norm of a vector
  template<typename T>
  inline T norm(const Vector<T> *v)
  {
    T result; nullify(result);
    for (T const* vIt = v->begin(); vIt != v->end(); ++vIt)
      result += *vIt * *vIt;
    return std::sqrt(result);
  }
  
  /// 2-norm of a vector
  template<typename T>
  inline T norm(const Vector<T>& v)
  {
    return norm(&v);
  }

  /// cross-product of two vectors: vector := vector x vector (only in 3d)
  template<typename T>
  void vectorProduct(const Vector<T>& x, 
		     const Vector<T>& y, 
		     Vector<T>& z)
  {
    FUNCNAME_DBG("vectorProduct()");
    TEST_EXIT_DBG(Global::getGeo(WORLD) == 3)("DIM_OF_WORLD != 3\n");
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
  }  

  // ---------------------------------------------------------------------------
  // Operations with WorldVector and WorldMatrix

  // times / divides
  // ---------------

  /// Scalar product: scalar := vector * vector
  template<typename T, typename S>
  inline typename traits::mult_type<T,S>::type 
  operator*(const WorldVector<T>& v1, const WorldVector<S>& v2) 
  {
    typename traits::mult_type<T,S>::type result;
    nullify(result);

    T const* v1It;
    S const* v2It;
    for (v1It = v1.begin(), v2It = v2.begin();
	 v1It != v1.end();
	 ++v1It, ++v2It)
      result += *v1It * *v2It;

    return result;
  }
  
  /// vector := vector * scalar (elementwise)
  template<typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, WorldVector<T> >::type
  operator*(WorldVector<T> const& v, S scal)
  {
    WorldVector<T> result = v;
    result *= scal; // calls operator*=(Vector<T>, S)
    return result;
  }

  /// vector := scalar * vector (elementwise)
  template<typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, WorldVector<T> >::type
  operator*(S scal, WorldVector<T> const& v)
  {
    WorldVector<T> result = v;
    result *= scal; // calls operator*=(Vector<T>, S)
    return result;
  }

  /// vector := vector / scalar (elementwise)
  template<typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, WorldVector<T> >::type
  operator/(WorldVector<T> const& v, S scal)
  {
    WorldVector<T> result = v;
    result /= scal;  // calls operator/=(Vector<T>, S)
    return result;
  }
  
  /// matrix *= scalar (elementwise)
//   template<typename T>
//   WorldMatrix<T>& operator*=(WorldMatrix<T>& m, T scal)
//   {
//     for (T* mIt = m.begin(); mIt != m.end(); mIt++)
//       *mIt *= scal;
// 
//     return m;
//   }
  
  /// matrix := matrix * scalar (elementwise)
  template <typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, WorldMatrix<T> >::type
  operator*(WorldMatrix<T> const& m, S scal)
  {
    WorldMatrix<T> result = m;
    result *= scal; // calls operator*=(Vector<T>, S)
    return result;
  }

  /// matrix := scalar * matrix (elementwise)
  template <typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, WorldMatrix<T> >::type
  operator*(S scal, WorldMatrix<T> const& m)
  {
    WorldMatrix<T> result = m;
    result *= scal; // calls operator*=(Vector<T>, S)
    return result;
  }

  /// matrix := matrix / scalar (elementwise)
  template <typename T, typename S>
  typename enable_if< traits::is_multiplicable<S, T>, WorldMatrix<T> >::type
  operator/(WorldMatrix<T> const& m, S scal)
  {
    WorldMatrix<T> result = m;
    result /= scal; // calls operator/=(Vector<T>, S)
    return result;
  }

  /// vector := matrix * vector
  template<typename T>
  WorldVector<T> operator*(const WorldMatrix<T>& M, const WorldVector<T>& v )
  {
    WorldVector<T> res;
    res.multMatrixVec(M,v);
    return res;
  }

  /// matrix := matrix * matrix
  template<typename T>
  WorldVector<WorldVector<T> > 
  operator*(const WorldVector<WorldVector<T> >& A, const WorldVector<WorldVector<T> >& B)
  {
    WorldVector<WorldVector<T> > result;
    nullify(result);
    for (size_t r = 0; r < num_rows(A); r++)
      for (size_t c = 0; c < num_cols(A); c++)
	for (size_t i = 0; i < num_cols(A); i++)
	  result[r][c] += A[r][i] * B[i][c];
    return result;
  }
  
  // plus / minus
  // ------------
  
// NOTE: call operators of Vector<T> directly
#if 0
  template<typename T>
  WorldVector<T>& operator+=(WorldVector<T>& v1,
			     const WorldVector<T>& v2)
  {
    add(v1, v2, v1);
    return v1;
  }

  template<typename T>
  WorldVector<T>& operator-=(WorldVector<T>& v1,
			     const WorldVector<T>& v2)
  {
    axpy(-1.0, v2, v1);
    return v1;
  }
#endif

  /// vector := vector + vector
  template <typename T, typename S>
  WorldVector<T>&
  operator+=(WorldVector<T>& v1, WorldVector<S> const& v2)
  {
    static_cast<Vector<T>&>(v1) += static_cast<Vector<S> const&>(v2);
    return v1;
  }
  
  /// vector := vector + vector
  template <typename T, typename S>
  WorldVector<typename traits::add_type<T, S>::type>
  operator+(WorldVector<T> result, const WorldVector<S>& v2)
  {
    result += v2; // calls operator+=(Vector<T>, Vector<T>)
    return result;
  }

  /// vector := vector - vector
  template <typename T, typename S>
  WorldVector<typename traits::add_type<T, S>::type>
  operator-(WorldVector<T> result, const WorldVector<S>& v2)
  {
    result -= v2; // calls operator-=(Vector<T>, Vector<T>)
    return result;
  }

  /// matrix += matrix
  template <typename T, typename S>
  WorldMatrix<T>& operator+=(WorldMatrix<T>& m1, const WorldMatrix<S>& m2)
  {
    T* m1It;
    S const* m2It;
    for (m1It = m1.begin(), m2It = m2.begin();
	 m1It != m1.end(); 
	 m1It++, m2It++)
      *m1It += *m2It;

    return m1;
  }

  /// matrix += matrix
  template <typename T, typename S>
  Matrix<T>& operator+=(Matrix<T>& m1, const Matrix<S>& m2)
  {
    T* m1It;
    S const* m2It;
    for (m1It = m1.begin(), m2It = m2.begin();
	 m1It != m1.end(); 
	 m1It++, m2It++)
      *m1It += *m2It;

    return m1;
  }

  /// matrix := matrix + matrix
  template <typename T, typename S>
  WorldMatrix<T> operator+(WorldMatrix<T> M1, const WorldMatrix<S>& M2 )
  {
    M1 += M2;
    return M1;
  }

  /// matrix := matrix + matrix
  template <typename T, typename S>
  Matrix<T> operator+(Matrix<T> M1, const Matrix<S>& M2 )
  {
    M1 += M2;
    return M1;
  }
  
  /// matrix -= matrix
  template <typename T, typename S>
  WorldMatrix<T>& operator-=(WorldMatrix<T>& m1, const WorldMatrix<S>& m2)
  {
    T *m1It;
    S const* m2It;
    for (m1It = m1.begin(), m2It = m2.begin();
	 m1It != m1.end(); 
	 m1It++, m2It++)
      *m1It -= *m2It;

    return m1;
  }
  
  /// matrix -= matrix
  template <typename T, typename S>
  Matrix<T>& operator-=(Matrix<T>& m1, const Matrix<S>& m2)
  {
    T *m1It;
    S const* m2It;
    for (m1It = m1.begin(), m2It = m2.begin();
	 m1It != m1.end(); 
	 m1It++, m2It++)
      *m1It -= *m2It;

    return m1;
  }

  /// matrix := matrix - matrix
  template <typename T, typename S>
  WorldMatrix<T> operator-(WorldMatrix<T> M1, const WorldMatrix<S>& M2 )
  {
    M1 -= M2;
    return M1;
  }

  /// matrix := matrix - matrix
  template <typename T, typename S>
  Matrix<T> operator-(Matrix<T> M1, const Matrix<S>& M2 )
  {
    M1 -= M2;
    return M1;
  }
  
  // unary minus operators
  // ---------------------

  /// vector := -vector (elementwise)
  template<typename T>
  WorldVector<T> operator-(WorldVector<T> v)
  {
    v *= -1.0;
    return v;
  }

  /// matrix := -matrix (elementwise)
  template<typename T>
  WorldMatrix<T> operator-(WorldMatrix<T> v)
  {
    v *= -1.0;
    return v;
  }
  
  // comparison operators
  // --------------------

  /// test for less-then (elementwise) up to DBL_TOL
  inline bool operator<(const WorldVector<double>& v1, const WorldVector<double>& v2) 
  {
    int dow = Global::getGeo(WORLD);
    for (int i = 0; i < dow; i++) {
      if (std::abs(v1[i] - v2[i]) < DBL_TOL) 
	continue;
      return v1[i] < v2[i];
    }
    return false;
  }

  /// test for equality (elementwise) up to DBL_TOL
  inline bool operator==(const WorldVector<double>& v1, const WorldVector<double>& v2) 
  {
    int dow = Global::getGeo(WORLD);
    for (int i = 0; i < dow; i++)
      if (std::abs(v1[i] - v2[i]) > DBL_TOL) 
	return false;

    return true;
  }
  

  /// Comparison operator for Vector<T>.
  template <typename T>
  inline bool operator==(Vector<T> const& lhs, Vector<T> const& rhs) 
  {
    if (lhs.getSize() != rhs.getSize()) 
      return false;

    T const* rhsIt;
    T const* lhsIt;
    for (rhsIt = rhs.begin(), lhsIt = lhs.begin();
	  rhsIt != rhs.end();
	  ++rhsIt, ++lhsIt)
      if (*lhsIt != *rhsIt) 
	return false;

    return true;
  }

  /// Comparison operator for Vector<T>.
  template <typename T>
  inline bool operator!=(Vector<T> const& lhs, Vector<T> const& rhs) 
  {
    return !(lhs == rhs);
  }
  
  

  /// Comparison operator for Matrix<T>.
  template <typename T>
  inline bool operator==(Matrix<T> const& lhs, Matrix<T> const& rhs) 
  {
    if (lhs.getNumRows() != rhs.getNumRows()) return false;
    if (lhs.getNumCols() != rhs.getNumCols()) return false;
    return (static_cast<Vector<T> const&>(lhs) == static_cast<Vector<T> const&>(rhs));
  }

  /// Comparison operator for Matrix<T>.
  template <typename T>
  inline bool operator!=(Matrix<T> const& lhs, Matrix<T> const& rhs) 
  {
    return !(lhs == rhs);
  }
  
  // special operators
  // -----------------
  
  /// wrapper for nullify
  template<typename T>
  void set_to_zero(WorldVector<WorldVector<T> >& mat)
  {
    nullify(mat);
  }

// NOTE: call norm(Vector<T>) directly
#if 0
  inline double norm(const WorldVector<double>& v)
  {
    double val = 0.0;
    for (int i = 0; i < Global::getGeo(WORLD); i++)
      val += v[i] * v[i];
    return sqrt(val);
  }
#endif

  /// returns the euclidian distance of a and b
  template<typename T, GeoIndex d>
  double absteukl(const FixVec<T,d>& a,const FixVec<T,d>& b)
  {
    double erg = 0.0;
    for (int i = 0; i < a.getSize(); ++i)
      erg += sqr(a[i] - b[i]);

    return std::sqrt(erg);
  }

} // end namespace AMDiS

#endif // AMDIS_MATVEC_OPERATIONS_H

