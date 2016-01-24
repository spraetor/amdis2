/** \file Math.h */

/** \defgroup Common Common
 */

#pragma once

#include <vector>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "AMDiS_fwd.h"
#include <traits/basic.hpp>
#include <traits/scalar_types.hpp>

namespace AMDiS
{
  namespace math
  {
    template <class T>
    Requires_t<concepts::Arithmetic<T>, T>
    constexpr abs(T a)
    {
      return  a >= 0 ? a : -a;
    }


    template <class T>
    Requires_t<concepts::Arithmetic<T>, T>
    constexpr sqr(T a)
    {
      return a*a;
    }


    template <class T0, class T1>
    Common_t<T0, T1>
    constexpr min(T0 a, T1 b)
    {
      return a > b ? b : a;
    }


    template <class T0, class T1>
    Common_t<T0, T1>
    constexpr max(T0 a, T1 b)
    {
      return a > b ? a : b;
    }

  } // end namespace math


  template <class T> inline void nullify(T& a)
  {
    a = 0;
  }

  template <class T> inline void nullify(std::vector<T>& a)
  {
    typename std::vector<T>::iterator it;
    for (it = a.begin(); it != a.end(); it++)
      nullify(*it);
  }

  template <class T> inline void nullify(mtl::dense_vector<T>& a)
  {
    T null;
    nullify(null);
    a = null;
  }

  /// Calculates factorial of i
  inline long fac(long i)
  {
    if (i <= 1)
      return 1;
    else
      return i * fac(i - 1);
  }

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
