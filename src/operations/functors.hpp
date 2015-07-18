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

#ifndef AMDIS_OPERATIONS_FUNCTOR_HPP
#define AMDIS_OPERATIONS_FUNCTOR_HPP

#include <complex>
#include <cmath>

#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/pow.hpp> 
#include "traits/types.hpp"
#include "operations/meta.hpp"

namespace AMDiS 
{
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
    template<typename T>
    struct identity : FunctorBase
    {
      typedef T result_type;
      int degree(int d0) const { return d0; }

      static T eval(const T& v) { return v; }
      T operator()(const T& v) const { return eval(v); }
    };

    /// constant(v) == val
    template<typename T>
    struct constant : FunctorBase
    {
      typedef T result_type;
      constant(T val_) : val(val_) {}

      template<typename V>
      result_type operator()(const V& v) const { return val; }

    private:
      T val;
    };
    
    template<typename T, typename S=T>
    struct add_constant : FunctorBase
    {
      typedef T result_type;
      S value;
      add_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v += value); }
    };
    
    template<typename T, typename S=T>
    struct minus_constant : FunctorBase
    {
      typedef T result_type;
      S value;
      minus_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v -= value); }
    };
    
    template<typename T, typename S=T>
    struct mult_constant : FunctorBase
    {
      typedef T result_type;
      S value;
      mult_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v *= value); }
    };
    
    template<typename T, typename S=T>
    struct div_constant : FunctorBase
    {
      typedef T result_type;
      S value;
      div_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v /= value); }
    };

    /// functor for operator+=
    template<typename T>
    struct assign : FunctorBase
    {
      typedef T result_type;
      
      static result_type& apply(T& v, T const& v0) { return (v = v0); }
      result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator+=
    template<typename T>
    struct add_assign : FunctorBase
    {
      typedef T result_type;
      
      static result_type& apply(T& v, T const& v0) { return (v += v0); }
      result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator*=
    template<typename T>
    struct mult_assign : FunctorBase
    {
      typedef T result_type;
      
      static result_type& apply(T& v, T const& v0) { return (v *= v0); }
      result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };
    

    /// abs(v) == |v|
    template<typename T>
    struct abs : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v) { return std::abs(v); }
      result_type operator()(const T &v) const { return eval(v); }
    };

    template<typename T>
    struct abs<std::complex<T> > : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v) { return std::norm(v); }
      result_type operator()(const T &v) const { return eval(v); }
    };

    /// max(a,b)
    template<typename T>
    struct max : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::max(v0, v1); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// min(a,b)
    template<typename T>
    struct min : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::min(v0, v1); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// max(|a|,|b|)
    template<typename T>
    struct abs_max : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::max(std::abs(v0), std::abs(v1)); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// min(|a|,|b|)
    template<typename T>
    struct abs_min : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::min(std::abs(v0), std::abs(v1)); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };
    

    /// apply a functor N times
    template<typename Functor, int N>
    struct apply
    {
      typedef typename Functor::result_type result_type;

      apply(const Functor& f_) : f(f_), inner(f_) {}
      int getDegree(int d0) const
      {
        return f.getDegree(inner.getDegree(d0));
      }

      template<typename V>
      static result_type eval(const V& v) { return Functor::eval(apply<Functor, N-1>::eval(v)); }
      template<typename V>
      result_type operator()(const V& v) const 
      {
	return f(inner(v));
      }

    private:
      const Functor& f;
      apply<Functor, N-1> inner; 
    };

    template<typename Functor>
    struct apply<Functor, 0>
    {
      typedef typename Functor::result_type result_type;

      apply(const Functor& f_) : f(f_) {}
      int getDegree(int d0) const { return d0; }

      template<typename V>
      static result_type eval(const V& v) { return v; }
      template<typename V>
      result_type operator()(const V& v) const { return v; }

    private:
      const Functor& f;
    };

    /// pow<p>(v) == v^p
    template<int p, typename T>
    struct pow : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return p*d0; }

      static result_type eval(const T& v) { return boost::math::pow<p>(v); }
      result_type operator()(const T& v) const { return eval(v); }
    };
  
    /// root<p>(v) == p-th-root(v)
    template<int p, typename T, typename Enabled = void>
    struct root_dispatch;

    template<int p, typename T>
    struct root : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return p*d0; } // optimal polynomial approximation degree ?

      static result_type eval(const T& v) { return root_dispatch<p,T>::eval(v); }
      result_type operator()(const T& v) const { return eval(v); }
    };

    template<int p, typename T, typename Enabled>
    struct root_dispatch { static T eval(const T& v) { return std::pow(v, 1.0/p); } };

    template<int p, typename T>
    struct root_dispatch<p, T, typename boost::enable_if<typename meta::is_power_of<p, 3>::type>::type> 
    {
      static T eval(const T& v) { return apply<root<3, T>, meta::log<p, 3>::value>::eval(v); }
    };
  
    template<int p, typename T>
    struct root_dispatch<p, T, typename boost::enable_if<typename meta::is_power_of<p, 2>::type>::type> 
    {
      static T eval(const T& v) { return apply<root<2, T>, meta::log<p, 2>::value>::eval(v); }
    };
  
    template<typename T>
    struct root_dispatch<3, T> { static T eval(const T& v) { return boost::math::cbrt(v); } };

    template<typename T>
    struct root_dispatch<2, T> { static T eval(const T& v) { return std::sqrt(v); } };
  
    template<typename T>
    struct root_dispatch<1, T> { static T eval(const T& v) { return v; } };
  
    template<typename T>
    struct root_dispatch<0, T> { static T eval(const T& v) { return 1.0; } };

  } // end namespace functors

} // end namespace AMDiS

#endif // AMDIS_OPERATIONS_FUNCTOR_HPP
