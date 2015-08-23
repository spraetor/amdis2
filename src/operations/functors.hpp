/** \file functors.hpp */

#pragma once

// std c++ headers
#include <complex>
#include <cmath>
#include <type_traits>

// boost headers
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/pow.hpp> 

// AMDiS headers
#include "traits/types.hpp"
#include "traits/scalar_types.hpp"
#include "traits/mult_type.hpp"
#include "operations/meta.hpp"

namespace AMDiS 
{
  struct FunctorBase
  {
    int getDegree() const { return 0; }
    
    template <class... Ints>
    typename enable_if< and_< std::is_integral<Ints>... >, int >::type
    getDegree(Ints...) const { return 0; }
  };
  
  
  template <class T, class Model>
  struct UnaryFunctor : FunctorBase
  {
    using result_type = decltype( Model::eval(std::declval<T>()) );
    constexpr result_type operator()(const T& v) const { return Model::eval(v); }
  };
  
  
  template <class T1, class T2, class Model>
  struct BinaryFunctor : FunctorBase
  {
    using result_type = decltype( Model::eval(std::declval<T1>(), std::declval<T2>()) );
    constexpr result_type operator()(const T1& v1, const T2& v2) const { return Model::eval(v1, v2); }
  };
    
  

  namespace functors
  {
    /// identity(v) == v
    template <class T>
    struct identity : FunctorBase
    {
      using result_type = T;
      constexpr int degree(int d0) const { return d0; }

      constexpr static T eval(const T& v) { return v; }
      constexpr T operator()(const T& v) const { return eval(v); }
    };

    /// constant(v) == val
    template <class T>
    struct constant : FunctorBase
    {
      using result_type = T;
      constant(T val_) : val(val_) {}

      template <class V>
      result_type operator()(V&&) const { return val; }

    private:
      T val;
    };
    
    /// ct_constant(v) == val
    template <class T, long val_>
    struct ct_constant : FunctorBase
    {
      using result_type = T;
      static constexpr T val = val_;
      
      template <class V> constexpr static T eval(V&&) { return val; }
      template <class V> constexpr static T apply(V&&) { return val; }
      template <class V> constexpr T operator()(V&&) const { return val; }
    };
    
    template <class T, class S=T>
    struct add_constant : FunctorBase
    {
      using result_type = T;
      S value;
      constexpr add_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v += value); }
    };
    
    template <class T, class S=T>
    struct minus_constant : FunctorBase
    {
      using result_type = T;
      S value;
      constexpr minus_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v -= value); }
    };
    
    template <class T, class S=T>
    struct mult_constant : FunctorBase
    {
      using result_type = T;
      S value;
      constexpr mult_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v *= value); }
    };
    
    template <class T, class S=T>
    struct div_constant : FunctorBase
    {
      using result_type = T;
      S value;
      constexpr div_constant(S value) : value(value) {}
      
      result_type& operator()(T& v) { return (v /= value); }
    };

    /// functor for operator+=
    template <class T>
    struct assign : FunctorBase
    {
      using result_type = T;
      
      constexpr static result_type& apply(T& v, T const& v0) { return (v = v0); }
      constexpr result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator+=
    template <class T>
    struct add_assign : FunctorBase
    {
      using result_type = T;
      
      constexpr static result_type& apply(T& v, T const& v0) { return (v += v0); }
      constexpr result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator*=
    template <class T>
    struct mult_assign : FunctorBase
    {
      using result_type = T;
      
      constexpr static result_type& apply(T& v, T const& v0) { return (v *= v0); }
      constexpr result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator/=
    template <class T>
    struct div_assign : FunctorBase
    {
      using result_type = T;
      
      constexpr static result_type& apply(T& v, T const& v0) { return (v /= v0); }
      constexpr result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for v=max(v,w)
    template <class T>
    struct max_assign : FunctorBase
    {
      using result_type = T;
      
      constexpr static result_type& apply(T& v, T const& v0) { return (v = std::max(v,v0)); }
      constexpr result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };

    /// functor for operator/=
    template <class T>
    struct min_assign : FunctorBase
    {
      using result_type = T;
      
      constexpr static result_type& apply(T& v, T const& v0) { return (v = std::min(v,v0)); }
      constexpr result_type& operator()(T& v, T const& v0) { return apply(v,v0); }
    };
    
    // -------------------------------------------------------------------------

    /// abs(v) == |v|
    template <class T>
    struct abs : FunctorBase
    {
      constexpr static int getDegree(int d0) { return d0; }
      constexpr static auto eval(const T& v) RETURNS( math::abs(v) )
      constexpr auto operator()(const T& v) const RETURNS( eval(v) )
    };

    template <class T>
    struct abs<std::complex<T>> : FunctorBase
    {
      constexpr static int getDegree(int d0) { return d0; }
      constexpr static auto eval(const T& v) RETURNS( std::norm(v) )
      constexpr auto operator()(const T& v) const RETURNS( eval(v) )
    };
    
    /// negate(v) == -v
    template <class T>
    struct negate : FunctorBase
    {
      constexpr static int getDegree(int d0) { return d0; }
      constexpr static auto eval(const T& v) RETURNS( -v )
      constexpr auto operator()(const T& v) const RETURNS( eval(v) )
    };

    /// a + b
    template <class T1, class T2 = T1>
    struct plus : FunctorBase
    {
      constexpr static int getDegree(int d0, int d1) { return math::max(d0, d1); }
      constexpr static auto eval(const T1& v0, const T2& v1) RETURNS( v0 + v1 )
      constexpr auto operator()(const T1& v1, const T2& v2) const RETURNS( eval(v1, v2) )
    };

    /// a - b
    template <class T1, class T2 = T1>
    struct minus : FunctorBase
    {
      constexpr static int getDegree(int d0, int d1) { return math::max(d0, d1); }
      constexpr static auto eval(const T1& v0, const T2& v1) RETURNS( v0 - v1 )
      constexpr auto operator()(const T1& v1, const T2& v2) const RETURNS( eval(v1, v2) )
    };

    /// a * b
    template <class T1, class T2 = T1>
    struct multiplies : FunctorBase
    {
      constexpr static int getDegree(int d0, int d1) { return d0 + d1; }
      constexpr static auto eval(const T1& v0, const T2& v1) RETURNS( v0 * v1 )
      constexpr auto operator()(const T1& v1, const T2& v2) const RETURNS( eval(v1, v2) )
    };

    /// a / b
    template <class T1, class T2 = T1>
    struct divides : FunctorBase
    {
      constexpr static int getDegree(int d0, int d1) { return d0 + d1; }
      constexpr static auto eval(const T1& v0, const T2& v1) RETURNS( v0 / v1 )
      constexpr auto operator()(const T1& v1, const T2& v2) const RETURNS( eval(v1, v2) )
    };

    /// max(a,b)
    template <class T>
    struct max : FunctorBase
    {
      constexpr static int getDegree(int d0, int d1) { return math::max(d0, d1); }
      constexpr static auto eval(const T& v0, const T& v1) RETURNS( math::max(v0, v1) )
      constexpr auto operator()(const T& v1, const T& v2) const RETURNS( eval(v1, v2) )
    };

    /// min(a,b)
    template <class T>
    struct min : FunctorBase
    {
      constexpr static int getDegree(int d0, int d1) { return math::max(d0, d1); }
      constexpr static auto eval(const T& v0, const T& v1) RETURNS( math::min(v0, v1) )
      constexpr auto operator()(const T& v1, const T& v2) const RETURNS( eval(v1, v2) )
    };

    /// max(|a|,|b|)
    template <class T>
    struct abs_max : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return d0; }

      constexpr static result_type eval(const T& v0, const T& v1) { return math::max(math::abs(v0), math::abs(v1)); }
      constexpr result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };

    /// min(|a|,|b|)
    template <class T>
    struct abs_min : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return d0; }

      constexpr static result_type eval(const T& v0, const T& v1) { return math::min(math::abs(v0), math::abs(v1)); }
      constexpr result_type operator()(const T &v0, const T& v1) const { return eval(v0,v1); }
    };    
    
    
    /// cross(v1, v2) = v1 x v2, TODO: find better name
    template <class T1, class T2>
    struct MyCross : FunctorBase
    {
      typedef typename traits::mult_type<T1, T2>::type value_type;
      constexpr int getDegree(int d, int d0, int d1) const { return d0+d1; }
      
      template <class Vec1, class Vec2>
      static value_type apply(size_t i, const Vec1 &v1, const Vec2& v2)
      {
        using size_type = Size_t<traits::category<Vec1>>;
        value_type result;
        
        TEST_EXIT_DBG( size(v1) == 3 && size(v1) == size(v2) )
          ("cross: inkompatible sizes!\n");

        size_type k = (i+1) % 3, l = (i+2) % 3;
        result = v1(k) * v2(l) - v1(l) * v2(k);
        return result;
      }
      
      template <class Vec1, class Vec2>
      value_type operator()(size_t i, const Vec1 &v1, const Vec2& v2) const 
      { 
        return apply(i, v1, v2); 
      }
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

    
    /// pow<p>(v) == v^p
    template <int p, class T>
    struct pow : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return p*d0; }

      constexpr static result_type eval(const T& v) { return boost::math::pow<p>(v); }
      constexpr result_type operator()(const T& v) const { return eval(v); }
    };
  
    /// root<p>(v) == p-th-root(v)
    template <int p, class T, class Enabled = void>
    struct root_dispatch;

    template <int p, class T>
    struct root : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return p*d0; } // optimal polynomial approximation degree ?

      constexpr static result_type eval(const T& v) { return root_dispatch<p,T>::eval(v); }
      constexpr result_type operator()(const T& v) const { return eval(v); }
    };

    template <int p, class T, class Enabled>
    struct root_dispatch { constexpr static T eval(const T& v) { return std::pow(v, 1.0/p); } };

    template <int p, class T>
    struct root_dispatch<p, T, typename enable_if<typename meta::is_power_of<p, 3>::type>::type> 
    {
      constexpr static T eval(const T& v) { return apply<root<3, T>, meta::log<p, 3>::value>::eval(v); }
    };
  
    template <int p, class T>
    struct root_dispatch<p, T, typename enable_if<typename meta::is_power_of<p, 2>::type>::type> 
    {
      constexpr static T eval(const T& v) { return apply<root<2, T>, meta::log<p, 2>::value>::eval(v); }
    };
  
    template <class T>
    struct root_dispatch<3, T> { constexpr static T eval(const T& v) { return boost::math::cbrt(v); } };

    template <class T>
    struct root_dispatch<2, T> { constexpr static T eval(const T& v) { return std::sqrt(v); } };
  
    template <class T>
    struct root_dispatch<1, T> { constexpr static T eval(const T& v) { return v; } };
  
    template <class T>
    struct root_dispatch<0, T> { constexpr static T eval(const T& v) { return 1.0; } };

  } // end namespace functors

} // end namespace AMDiS
