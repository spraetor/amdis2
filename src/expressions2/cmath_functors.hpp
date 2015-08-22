/** \file cmath_functors.hpp */

#pragma once

#include <cmath>
#include <type_traits>

#include <Math.h>
#include <operations/functors.hpp>

namespace AMDiS
{
  namespace functors
  {
    template <class T>
    struct Signum : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 3*d0; }

      constexpr static result_type eval(const T& v) { return (v > T(0) ? T(1) : T(-1)); }
      constexpr result_type operator()(const T& v) const { return eval(v); }
    };
        
    
    /// Expressions for ceiling a value
    template <class T>
    struct Ceil : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return d0; }

      constexpr static result_type eval(const T& v) { return std::ceil(v); }
      constexpr result_type operator()(const T& v) const { return eval(v); }
    };
    
    
    /// Expressions to round to a lower integer
    template <class T>
    struct Floor : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return d0; }

      constexpr static result_type eval(const T& v) { return std::floor(v); }
      constexpr result_type operator()(const T& v) const { return eval(v); }
    };
    
    
    // -------------------------------------------------------------------------
    
    
    /// Expressions for the exponential function
    template <class T>
    struct Exp : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }
      
      constexpr static result_type eval(const T& v) { return std::exp(v); }
      constexpr result_type operator()(const T& v) const { return eval(v); }
    };
    

    /// Expression for the logarithm
    template <class T>
    struct Log : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }
      
      constexpr static result_type eval(const T& v) { return std::log(v); }
      constexpr result_type operator()(const T& v) const { return eval(v); }
    };
    
    
    // -------------------------------------------------------------------------
    
    
    /// sin(v)
    template <class T>
    struct Sin : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::sin(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    /// cos(v)
    template <class T>
    struct Cos : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::cos(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    /// tan(v)
    template <class T>
    struct Tan : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::tan(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    // -------------------------------------------------------------------------
    
    
    /// asin(v)
    template <class T>
    struct Asin : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::asin(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    /// acos(v)
    template <class T>
    struct Acos : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::acos(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    /// atan(v)
    template <class T>
    struct Atan : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::atan(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    /// atan(v)
    template <class T1, class T2>
    struct Atan2 : FunctorBase
    {
      using result_type = typename std::common_type<T1, T2>::type;
      constexpr int getDegree(int d0, int d1) const { return 2*(d0+d1); }

      constexpr static result_type eval(const T1& v1, const T2& v2) { return std::atan2(v1, v2); }
      constexpr result_type operator()(const T1& v1, const T2& v2) const { return eval(v1, v2); }
    };
    
    
    // -------------------------------------------------------------------------
    
    
    /// sinh(v)
    template <class T>
    struct Sinh : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::sinh(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    /// cosh(v)
    template <class T>
    struct Cosh : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::cosh(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    /// tanh(v)
    template <class T>
    struct Tanh : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::tanh(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    // -------------------------------------------------------------------------
    
    
    /// sinh(v)
    template <class T>
    struct Asinh : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::asinh(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    /// cosh(v)
    template <class T>
    struct Acosh : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::acosh(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
    
    /// tanh(v)
    template <class T>
    struct Atanh : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const { return 2*d0; }

      constexpr static result_type eval(const T& v) { return std::atanh(v); }
      constexpr result_type operator()(const T &v) const { return eval(v); }
    };
    
  } // end namespace functors
} // end namespace AMDiS
