/** \file functors.hpp */

#pragma once

// std c++ headers
#include <cmath>
#include <complex>
#include <type_traits>

// boost headers
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/pow.hpp>

// AMDiS headers
#include "operations/functor_generator.hpp"
#include "operations/meta.hpp"
#include "traits/basic.hpp"
#include "traits/scalar_types.hpp"

namespace AMDiS
{
  namespace functors
  {
    /// identity(v) == v
    AMDIS_MAKE_UNARY_FUNCTOR( identity , d0 , v )

    /// constant(v) == val
    template <class T>
    struct constant : FunctorBase
    {
      constant(T val_) : val(val_) {}

      template <class V>
      T operator()(V&&) const
      {
        return val;
      }

    private:
      T val;
    };

    /// ct_constant(v) == val
    template <class T, long val_>
    struct ct_constant : FunctorBase
    {
      static constexpr T val = val_;

      template <class V> static constexpr T eval(V&&)
      {
        return val;
      }
      template <class V> static constexpr T apply(V&&)
      {
        return val;
      }
      template <class V> constexpr T operator()(V&&) const
      {
        return val;
      }
    };

    /// abs(v) == |v|
    template <class T>
    struct abs : FunctorBase
    {
      static constexpr int getDegree(int d0)
      {
        return d0;
      }
      static constexpr auto eval(const T& v) RETURNS( math::abs(v) )
      constexpr auto operator()(const T& v) const RETURNS( eval(v) )
    };

    // specialization of abs for complex values
    template <class T>
    struct abs<std::complex<T>> : FunctorBase
    {
      static constexpr int getDegree(int d0)
      {
        return d0;
      }
      static constexpr auto eval(const T& v) RETURNS( std::norm(v) )
      constexpr auto operator()(const T& v) const RETURNS( eval(v) )
    };

    /// negate(v) == -v
    AMDIS_MAKE_UNARY_FUNCTOR( negate, d0, -v  )

    AMDIS_MAKE_BINARY_FUNCTOR( plus,  math::max(d0, d1), v0 + v1  )
    AMDIS_MAKE_BINARY_FUNCTOR( minus, math::max(d0, d1), v0 - v1  )
    AMDIS_MAKE_BINARY_FUNCTOR( multiplies, d0 + d1,      v0 * v1  )
    AMDIS_MAKE_BINARY_FUNCTOR( divides,    d0 + d1,      v0 / v1  )


    // _____ logical functors _________________________________________________
    
    AMDIS_MAKE_BINARY_FUNCTOR( equal,   0, v0 == v1  )
    AMDIS_MAKE_BINARY_FUNCTOR( unequal, 0, v0 != v1  )
    AMDIS_MAKE_BINARY_FUNCTOR( less,    0, v0 < v1  )
    AMDIS_MAKE_BINARY_FUNCTOR( greater, 0, v0 > v1  )
    
    AMDIS_MAKE_BINARY_FUNCTOR( logical_and, 0, v0 && v1  )
    AMDIS_MAKE_BINARY_FUNCTOR( logical_or,  0, v0 || v1  )
    
    AMDIS_MAKE_BINARY_FUNCTOR( max, math::max(d0, d1), math::max(v0, v1)  )
    AMDIS_MAKE_BINARY_FUNCTOR( min, math::max(d0, d1), math::min(v0, v1)  )


    /// max(|a|,|b|)
    AMDIS_MAKE_BINARY_FUNCTOR( abs_max, math::max(d0, d1), math::max(math::abs(v0), math::abs(v1))  )

    /// min(|a|,|b|)
    AMDIS_MAKE_BINARY_FUNCTOR( abs_min, math::max(d0, d1), math::min(math::abs(v0), math::abs(v1))  )


    /// conditional(a,b,c) = a ? b : c
    template <class T1, class T2>
    struct conditional : FunctorBase
    {
      constexpr int getDegree(int d0, int d1, int d2) const
      {
        return math::max(d1, d2);
      }
      static constexpr typename std::common_type<T1, T2>::type 
      eval(bool cond, T1 const& v1, T2 const& v2)
      {
        return cond ? v1 : v2;
      }
      constexpr auto operator()(bool cond, T1 const& v1, T2 const& v2) const RETURNS
      ( 
        eval(cond, v1, v2) 
      )
    };


    /// cross(v1, v2) = v1 x v2, TODO: find better name
    template <class T1, class T2>
    struct MyCross : FunctorBase
    {
      using value_type = decltype( std::declval<T1>() * std::declval<T2>() );
      constexpr int getDegree(int /*d*/, int d0, int d1) const
      {
        return d0+d1;
      }

      template <class Vec1, class Vec2>
      static value_type eval(size_t i, const Vec1& v1, const Vec2& v2)
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
      auto operator()(size_t i, const Vec1& v1, const Vec2& v2) const RETURNS
      (
        eval(i, v1, v2)
      )
    };


    /// apply a functor N times
    template <class Functor, int N>
    struct apply
    {
      apply(Functor const& f_) : f(f_), inner(f_) {}

      int getDegree(int d0) const
      {
        return f.getDegree(inner.getDegree(d0));
      }

      template <class V>
      static auto eval(const V& v) RETURNS
      (
        Functor::eval(apply<Functor, N-1>::eval(v))
      )

      template <class V>
      auto operator()(const V& v) const RETURNS( f(inner(v)) )

    private:
      Functor f;
      apply<Functor, N-1> inner;
    };

    template <class Functor>
    struct apply<Functor, 0>
    {
      apply(Functor const& f_) : f(f_) {}
      int getDegree(int d0) const
      {
        return d0;
      }

      template <class V>
      static auto eval(const V& v) RETURNS( v )
      template <class V>
      auto operator()(const V& v) const RETURNS( v )

    private:
      Functor f;
    };



    // -------------------------------------------------------------------------

    template <class F, int arg, class G>
    struct compose;

    template <class F, class G>
    struct compose<F, 1, G>
    {
      template <class T>
      auto operator()(T const& v, T const& v0) RETURNS( f(g(v), v0) )
      
    private:
      F f;
      G g;
    };

    template <class F, class G>
    struct compose<F, 2, G>
    {
      template <class T>
      auto operator()(T const& v, T const& v0) RETURNS( f(v, g(v0)) )
      
    private:
      F f;
      G g;
    };


    /// pow<p>(v) == v^p
    template <int p, class T>
    struct pow : FunctorBase
    {
      constexpr int getDegree(int d0) const
      {
        return p*d0;
      }

      static constexpr T eval(const T& v)
      {
        return boost::math::pow<p>(v);
      }
      constexpr T operator()(const T& v) const
      {
        return eval(v);
      }
    };

    /// root<p>(v) == p-th-root(v)
    template <int p, class T, class = void>
    struct root_dispatch;

    template <int p, class T>
    struct root : FunctorBase
    {
      constexpr int getDegree(int d0) const
      {
        return p*d0;    // optimal polynomial approximation degree ?
      }

      static constexpr T eval(const T& v)
      {
        return root_dispatch<p,T>::eval(v);
      }
      constexpr T operator()(const T& v) const
      {
        return eval(v);
      }
    };

    template <int p, class T, class>
    struct root_dispatch
    {
      static constexpr T eval(const T& v)
      {
        return std::pow(v, 1.0/p);
      }
    };

    template <int p, class T>
    struct root_dispatch<p, T,
      Requires_t<meta::is_power_of<p, 3>> >
    {
      static constexpr T eval(const T& v)
      {
        return apply<root<3, T>, meta::log<p, 3>::value>::eval(v);
      }
    };

    template <int p, class T>
    struct root_dispatch<p, T,
      Requires_t<meta::is_power_of<p, 2>> >
    {
      static constexpr T eval(const T& v)
      {
        return apply<root<2, T>, meta::log<p, 2>::value>::eval(v);
      }
    };

    template <class T>
    struct root_dispatch<3, T>
    {
      static constexpr T eval(const T& v)
      {
        return boost::math::cbrt(v);
      }
    };

    template <class T>
    struct root_dispatch<2, T>
    {
      static constexpr T eval(const T& v)
      {
        return std::sqrt(v);
      }
    };

    template <class T>
    struct root_dispatch<1, T>
    {
      static constexpr T eval(const T& v)
      {
        return v;
      }
    };

    template <class T>
    struct root_dispatch<0, T>
    {
      static constexpr T eval(const T& /*v*/)
      {
        return 1.0;
      }
    };

  } // end namespace functors

} // end namespace AMDiS
