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



/** \file functors.hpp */

#pragma once

#include <complex>
#include <cmath>

#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/pow.hpp> 

#include "operations/meta.hpp"
#include "traits/mult_type.hpp"

namespace AMDiS 
{
  /// base class for all functors 
  struct FunctorBase
  {
    int getDegree() const { return 0; }
    int getDegree(int d0) const { return 0; }
    int getDegree(int d0, int d1) const { return 0; }
    int getDegree(int d0, int d1, int d2) const { return 0; }
    int getDegree(int d0, int d1, int d2, int d3) const { return 0; }
  };

  namespace functors
  {
    /// identity(v) == v
    template <class T>
    struct identity : FunctorBase
    {
      typedef T result_type;
      int degree(int d0) const { return d0; }

      static T eval(const T& v) { return v; }
      static T apply(const T& v) { return v; }
      T operator()(const T& v) const { return eval(v); }
    };

    // -------------------------------------------------------------------------
    /// constant(v) == val
    template <class T>
    struct constant : FunctorBase
    {
      typedef T result_type;
      constant(T val_) : val(val_) {}

      template <class V>
      result_type operator()(const V& v) const { return val; }

    private:
      T val;
    };
    
    /// ct_constant(v) == val
    template <class T, long val_>
    struct ct_constant : FunctorBase
    {
      typedef T result_type;
      static constexpr T val = val_;
      
      template <class V>
      static T eval(const V& v) { return val; }
      template <class V>
      static T apply(const V& v) { return val; }
      template <class V>
      T operator()(const V& v) const { return val; }
    };
    

    // -------------------------------------------------------------------------
    /// abs(v) == |v|
    template <class T>
    struct abs : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v) { return std::abs(v); }
      static result_type apply(const T& v) { return std::abs(v); }
      result_type operator()(const T &v) const { return eval(v); }
    };

    /// \cond HIDDEN_SYMBOLS
    template <class T>
    struct abs<std::complex<T> > : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v) { return std::norm(v); }
      static result_type apply(const T& v) { return std::norm(v); }
      result_type operator()(const T &v) const { return eval(v); }
    };
    /// \endcond

    /// max(a,b)
    template <class T>
    struct max : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::max(v0, v1); }
      static result_type apply(T& v0, const T& v1) { return (v0 = std::max(v0, v1)); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// min(a,b)
    template <class T>
    struct min : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::min(v0, v1); }
      static result_type apply(T& v0, const T& v1) { return (v0 = std::min(v0, v1)); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// max(|a|,|b|)
    template <class T>
    struct abs_max : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::max(std::abs(v0), std::abs(v1)); }
      static result_type apply(T& v0, const T& v1) { return (v0 = std::max(std::abs(v0), std::abs(v1))); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// min(|a|,|b|)
    template <class T>
    struct abs_min : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::min(std::abs(v0), std::abs(v1)); }
      static result_type apply(T& v0, const T& v1) { return (v0 = std::min(std::abs(v0), std::abs(v1))); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };
    
    
    /// cross(v1, v2) = v1 x v2
    template <class T1, class T2>
    struct MyCross : FunctorBase
    {
      typedef typename traits::mult_type<T1, T2>::type value_type;
      int getDegree(int d, int d0, int d1) const { return d0+d1; }
      
      template <class Vec1, class Vec2>
      static value_type apply(size_t i, const Vec1 &v1, const Vec2& v2)
      {
	typedef typename traits::category<Vec1>::size_type size_type;
	value_type result;
	
	assert( size(v1) == 3 && size(v1) == size(v2) ); // ("cross: inkompatible sizes!");

	size_type k = (i+1) % 3, l = (i+2) % 3;
	result = v1(k) * v2(l) - v1(l) * v2(k);
	return result;
      }
      
      template <class Vec1, class Vec2>
      value_type operator()(size_t i, const Vec1 &v1, const Vec2& v2) const { return apply(i, v1, v2); }
    };
    

    // -------------------------------------------------------------------------
    /// apply a functor \p N times
    template <class Functor, int N>
    struct apply
    {
      typedef typename Functor::result_type result_type;

      apply(const Functor& f_) : f(f_), inner(f_) {}
      int getDegree(int d0) const
      {
        return f.getDegree(inner.getDegree(d0));
      }

      template <class V>
      static result_type eval(const V& v) { return Functor::eval(apply<Functor, N-1>::eval(v)); }
      template <class V>
      result_type operator()(const V& v) const 
      {
	return f(inner(v));
      }

    private:
      const Functor& f;
      apply<Functor, N-1> inner; 
    };

    /// \cond HIDDEN_SYMBOLS
    template <class Functor>
    struct apply<Functor, 0>
    {
      typedef typename Functor::result_type result_type;

      apply(const Functor& f_) : f(f_) {}
      int getDegree(int d0) const { return d0; }

      template <class V>
      static result_type eval(const V& v) { return v; }
      template <class V>
      result_type operator()(const V& v) const { return v; }

    private:
      const Functor& f;
    };
    /// \endcond
    

    // -------------------------------------------------------------------------
    /// pow<p>(v) == v^p
    template <int p, class T>
    struct pow : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return p*d0; }

      static result_type eval(const T& v) { return boost::math::pow<p>(v); }
      result_type operator()(const T& v) const { return eval(v); }
    };
  
    // -------------------------------------------------------------------------
    /// \cond HIDDEN_SYMBOLS
    template <int p, class T, class Enabled = void>
    struct root_dispatch;
    /// \endcond

    /// root<p>(v) == p-th-root(v)
    template <int p, class T>
    struct root : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return p*d0; } // optimal polynomial approximation degree ?

      static result_type eval(const T& v) { return root_dispatch<p,T>::eval(v); }
      static result_type apply(const T& v) { return root_dispatch<p,T>::eval(v); }
      result_type operator()(const T& v) const { return eval(v); }
    };
    
    /// sqrt(v) == square-root-of(v)
    template <class T>
    struct sqrt : root<2,T> {};
    

    /// \cond HIDDEN_SYMBOLS
    template <int p, class T, class Enabled>
    struct root_dispatch { static T eval(const T& v) { return std::pow(v, 1.0/p); } };

    template <int p, class T>
    struct root_dispatch<p, T, typename enable_if<meta::is_power_of<p, 3> >::type > 
    {
      static T eval(const T& v) { return apply<root<3, T>, meta::log<p, 3>::value>::eval(v); }
    };
  
    template <int p, class T>
    struct root_dispatch<p, T, typename enable_if<meta::is_power_of<p, 2> >::type > 
    {
      static T eval(const T& v) { return apply<root<2, T>, meta::log<p, 2>::value>::eval(v); }
    };
  
    template <class T>
    struct root_dispatch<3, T> { static T eval(const T& v) { return boost::math::cbrt(v); } };

    template <class T>
    struct root_dispatch<2, T> { static T eval(const T& v) { return std::sqrt(v); } };
  
    template <class T>
    struct root_dispatch<1, T> { static T eval(const T& v) { return v; } };
  
    template <class T>
    struct root_dispatch<0, T> { static T eval(const T& v) { return 1.0; } };
    /// \endcond
    
    
    
    // -------------------------------------------------------------------------
    
    template <class F, int arg, class G>
    struct compose;
    
    template <class F, class G>
    struct compose<F, 1, G>
    {
      typedef typename F::result_type result_type;
      F f;
      G g;
      
      template <class T>
      result_type& operator()(T& v, T const& v0) { return f(g(v), v0); }
    };
    
    template <class F, class G>
    struct compose<F, 2, G>
    {
      typedef typename F::result_type result_type;
      F f;
      G g;
      
      template <class T>
      result_type& operator()(T& v, T const& v0) { return f(v, g(v0)); }
    };

  } // end namespace functors

} // end namespace AMDiS
