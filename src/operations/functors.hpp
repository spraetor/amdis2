/** \file functors.hpp */

#pragma once

#include <complex>
#include <cmath>
#include <type_traits>

#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/pow.hpp> 

#include "traits/types.hpp"
#include "traits/scalar_types.hpp"
#include "operations/meta.hpp"

namespace AMDiS 
{
  struct FunctorBase
  {
    int getDegree() const { return 0; }
    
    template <class... Ints>
    typename enable_if< and_< std::is_integral<Ints>... >, int >::type
    getDegree(Ints&&...) const { return 0; }
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
      T operator()(const T& v) const { return eval(v); }
    };

    /// constant(v) == val
    template <class T>
    struct constant : FunctorBase
    {
      typedef T result_type;
      constant(T val_) : val(val_) {}

      template <class V>
      result_type operator()(V&&) const { return val; }

    private:
      T val;
    };
    
    template <class T, class S=T>
    struct add_constant : FunctorBase
    {
      typedef T result_type;
      S value;
      add_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v += value); }
    };
    
    template <class T, class S=T>
    struct minus_constant : FunctorBase
    {
      typedef T result_type;
      S value;
      minus_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v -= value); }
    };
    
    template <class T, class S=T>
    struct mult_constant : FunctorBase
    {
      typedef T result_type;
      S value;
      mult_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v *= value); }
    };
    
    template <class T, class S=T>
    struct div_constant : FunctorBase
    {
      typedef T result_type;
      S value;
      div_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v /= value); }
    };

    /// functor for operator+=
    template <class T>
    struct assign : FunctorBase
    {
      typedef T result_type;
      
      static result_type& apply(T& v, T const& v0) { return (v = v0); }
      result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator+=
    template <class T>
    struct add_assign : FunctorBase
    {
      typedef T result_type;
      
      static result_type& apply(T& v, T const& v0) { return (v += v0); }
      result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator*=
    template <class T>
    struct mult_assign : FunctorBase
    {
      typedef T result_type;
      
      static result_type& apply(T& v, T const& v0) { return (v *= v0); }
      result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator/=
    template <class T>
    struct div_assign : FunctorBase
    {
      typedef T result_type;
      
      static result_type& apply(T& v, T const& v0) { return (v /= v0); }
      result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for v=max(v,w)
    template <class T>
    struct max_assign : FunctorBase
    {
      typedef T result_type;
      
      static result_type& apply(T& v, T const& v0) { return (v = std::max(v,v0)); }
      result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator/=
    template <class T>
    struct min_assign : FunctorBase
    {
      typedef T result_type;
      
      static result_type& apply(T& v, T const& v0) { return (v = std::min(v,v0)); }
      result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };
    

    /// abs(v) == |v|
    template <class T>
    struct abs : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v) { return std::abs(v); }
      result_type operator()(const T &v) const { return eval(v); }
    };

    template <class T>
    struct abs<std::complex<T> > : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v) { return std::norm(v); }
      result_type operator()(const T &v) const { return eval(v); }
    };

    /// max(a,b)
    template <class T>
    struct max : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::max(v0, v1); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// min(a,b)
    template <class T>
    struct min : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::min(v0, v1); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// max(|a|,|b|)
    template <class T>
    struct abs_max : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::max(std::abs(v0), std::abs(v1)); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// min(|a|,|b|)
    template <class T>
    struct abs_min : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return d0; }

      static result_type eval(const T& v0, const T& v1) { return std::min(std::abs(v0), std::abs(v1)); }
      result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };
    

    /// apply a functor N times
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
      result_type operator()(const V& v) const { return f(inner(v)); }

    private:
      const Functor& f;
      apply<Functor, N-1> inner; 
    };

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

    /// pow<p>(v) == v^p
    template <int p, class T>
    struct pow : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return p*d0; }

      static result_type eval(const T& v) { return boost::math::pow<p>(v); }
      result_type operator()(const T& v) const { return eval(v); }
    };
  
    /// root<p>(v) == p-th-root(v)
    template <int p, class T, class Enabled = void>
    struct root_dispatch;

    template <int p, class T>
    struct root : FunctorBase
    {
      typedef T result_type;
      int getDegree(int d0) const { return p*d0; } // optimal polynomial approximation degree ?

      static result_type eval(const T& v) { return root_dispatch<p,T>::eval(v); }
      result_type operator()(const T& v) const { return eval(v); }
    };

    template <int p, class T, class Enabled>
    struct root_dispatch { static T eval(const T& v) { return std::pow(v, 1.0/p); } };

    template <int p, class T>
    struct root_dispatch<p, T, typename enable_if<typename meta::is_power_of<p, 3>::type>::type> 
    {
      static T eval(const T& v) { return apply<root<3, T>, meta::log<p, 3>::value>::eval(v); }
    };
  
    template <int p, class T>
    struct root_dispatch<p, T, typename enable_if<typename meta::is_power_of<p, 2>::type>::type> 
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

  } // end namespace functors

} // end namespace AMDiS
