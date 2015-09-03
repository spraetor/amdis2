/** \file DOFVectorOperations.hpp */

#pragma once

#include "DOFVector.h"

namespace AMDiS
{
  template <class T>
  inline double norm(DOFVector<T> const& vec)
  {
    return vec.nrm2();
  }

  template <class T>
  inline double L2Norm(DOFVector<T> const& vec)
  {
    return vec.L2Norm();
  }

  template <class T>
  inline double H1Norm(DOFVector<T> const& vec)
  {
    return vec.H1Norm();
  }

  template <class T>
  inline void print(DOFVector<T> const& vec)
  {
    vec.print();
  }

  template <class T>
  inline void set_to_zero(DOFVector<T>& v)
  {
    T my_zero;
    nullify(my_zero);
    std::fill(v.begin(), v.end(), my_zero);
  }

  inline double min(DOFVector<double> const& v)
  {
    return v.min();
  }

  inline double max(DOFVector<double> const& v)
  {
    return v.max();
  }
  
  template <class T>
  T integrate(DOFVector<T> const& vec)
  {
    return vec.Int();
  }


  /* ----- OPERATORS WITH DOFVECTORS ---------------------------------------- */


  template <class T>
  DOFVector<T> operator*(DOFVector<T> v, T d)
  {
    v *= d;
    return v;
  }

  template <class T>
  DOFVector<T> operator*(T d, DOFVector<T> v)
  {
    v *= d;
    return v;
  }


  // scalar product
  template <class T>
  T operator*(DOFVector<T> const& x, DOFVector<T> const& y)
  {
    FUNCNAME("operator*(DOFVector<T>& x, DOFVector<T>& y)");
    const DOFAdmin* admin = NULL;

    TEST_EXIT(x.getFeSpace() && y.getFeSpace())
    ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT((admin = x.getFeSpace()->getAdmin()) && (admin == y.getFeSpace()->getAdmin()))
    ("no admin or different admins: %8X, %8X\n",
     x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT(x.getSize() == y.getSize())("different sizes\n");

    T dot;
    nullify(dot);

    DOFConstIterator<T> xIterator(&x, USED_DOFS);
    DOFConstIterator<T> yIterator(&y, USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end(); ++xIterator, ++yIterator)
      dot += (*xIterator) * (*yIterator);

    return dot;
  }


  // addition
  template <class T>
  DOFVector<T> operator+(DOFVector<T> v1 , const DOFVector<T>& v2)
  {
    v1 += v2;
    return v1;
  }

  // subtraction
  template <class T>
  DOFVector<T> operator-(DOFVector<T> v1 , const DOFVector<T>& v2)
  {
    v1 -= v2;
    return v1;
  }

} // end namespace AMDiS
