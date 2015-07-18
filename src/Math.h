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


/** \file Math.h */


/** \defgroup Common Common
 */

#ifndef AMDIS_MATH_H
#define AMDIS_MATH_H

#include <vector>
#include "AMDiS_fwd.h"
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace AMDiS {
  
  // ===== some simple template functions ====================================== 
#if 0
  template <class T> inline T abs(T a) 
  {
    return (a >= 0 ? a : -a);
  }
#endif

  template <class T> inline T sqr(T a) 
  {
    return a*a;
  }

  template <class T> inline void nullify(T &a)
  {
    a = 0;
  }

  template <class T> inline void nullify(std::vector<T> &a)
  {
    typename std::vector<T>::iterator it;
    for (it = a.begin(); it != a.end(); it++)
      nullify(*it);
  }

  template <class T> inline void nullify(mtl::dense_vector<T> &a)
  {
    T null;
    nullify(null);
    a = null;
  }

  template <class T> inline void nullify(WorldVector<T> &a)
  {
    T null; nullify(null);
    a = null;
  }

  template <class T> inline void nullify(WorldMatrix<T> &a)
  {
    T null; nullify(null);
    a = null;
  }

  /// Calculates factorial of i
  int fac(int i);
  
  /// check for inf and nan values
  inline bool isNumber(double val)
  {
    return !boost::math::isnan(val) && !boost::math::isinf(val);
  }

  // ===== some predefined values ===============================================
  constexpr double m_e = 	2.7182818284590452354;
  constexpr double m_log2e = 	1.4426950408889634074;
  constexpr double m_log10e = 	0.43429448190325182765;
  constexpr double m_ln2 = 	0.69314718055994530942;
  constexpr double m_ln10 = 	2.30258509299404568402;
  constexpr double m_pi = 	3.14159265358979323846;
  constexpr double m_pi_2 = 	1.57079632679489661923;
  constexpr double m_pi_4 = 	0.78539816339744830962;
  constexpr double m_1_pi = 	0.31830988618379067154;
  constexpr double m_2_pi = 	0.63661977236758134308;
  constexpr double m_2_sqrtpi = 1.12837916709551257390;
  constexpr double m_sqrt2 = 	1.41421356237309504880;
  constexpr double m_sqrt1_2 = 	0.70710678118654752440;

  // ===== tolerance for floating point comparison ==============================
#define DBL_TOL DBL_EPSILON
#define FLT_TOL FLT_EPSILON


} // end namespace AMDiS