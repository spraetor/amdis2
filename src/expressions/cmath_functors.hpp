#pragma once

// std c++ headers
#include <cmath>
#include <type_traits>

// AMDiS includes
#include "Math.hpp"
#include "operations/functors.hpp"

/// Macro that generates a unary functor.
/**
 *  \p NAME    Name of the class.
 *  \p DEGREE  Expression in 'd0' that gives the polynomial degree, with
 *             'd0' the degree of the argument passed to the functor.
 *  \p FCT     Name of a unary c++-function that represents the functor.
 */
#define AMDIS_MAKE_UNARY_FUNCTOR( NAME, DEGREE, FCT )     \
    template <class T>                                    \
    struct NAME : FunctorBase                             \
    {                                                     \
      using result_type = T;                              \
      constexpr int getDegree(int d0) const               \
      {                                                   \
        return DEGREE ;                                   \
      }                                                   \
      static constexpr result_type eval(T const& v)       \
      {                                                   \
        return FCT(v) ;                                   \
      }                                                   \
      constexpr result_type operator()(T const& v) const  \
      {                                                   \
        return eval(v);                                   \
      }                                                   \
    };


namespace AMDiS
{
  namespace functors
  {
    /// Functor that represents the signum function
    template <class T>
    struct Signum : FunctorBase
    {
      using result_type = T;
      constexpr int getDegree(int d0) const
      {
        return 3*d0;
      }

      static constexpr result_type eval(T const& v)
      {
        return (v > T{0} ? T{1} : T{-1});
      }
      constexpr result_type operator()(T const& v) const
      {
        return eval(v);
      }
    };

    // generated unary functors using a macro...

    AMDIS_MAKE_UNARY_FUNCTOR( Ceil,    d0, std::ceil  )
    AMDIS_MAKE_UNARY_FUNCTOR( Floor,   d0, std::floor )
    AMDIS_MAKE_UNARY_FUNCTOR( Exp,   2*d0, std::exp   )
    AMDIS_MAKE_UNARY_FUNCTOR( Log,   2*d0, std::log   )
    AMDIS_MAKE_UNARY_FUNCTOR( Sin,   2*d0, std::sin   )
    AMDIS_MAKE_UNARY_FUNCTOR( Cos,   2*d0, std::cos   )
    AMDIS_MAKE_UNARY_FUNCTOR( Tan,   2*d0, std::tan   )
    AMDIS_MAKE_UNARY_FUNCTOR( Asin,  2*d0, std::asin  )
    AMDIS_MAKE_UNARY_FUNCTOR( Acos,  2*d0, std::acos  )
    AMDIS_MAKE_UNARY_FUNCTOR( Atan,  2*d0, std::atan  )
    AMDIS_MAKE_UNARY_FUNCTOR( Sinh,  2*d0, std::sinh  )
    AMDIS_MAKE_UNARY_FUNCTOR( Cosh,  2*d0, std::cosh  )
    AMDIS_MAKE_UNARY_FUNCTOR( Tanh,  2*d0, std::tanh  )
    AMDIS_MAKE_UNARY_FUNCTOR( Asinh, 2*d0, std::asinh )
    AMDIS_MAKE_UNARY_FUNCTOR( Acosh, 2*d0, std::acosh )
    AMDIS_MAKE_UNARY_FUNCTOR( Atanh, 2*d0, std::atanh )


    /// Binary functor representing the cmath function std::atan2(v)
    template <class T1, class T2>
    struct Atan2 : FunctorBase
    {
      using result_type = typename std::common_type<T1, T2>::type;
      constexpr int getDegree(int d0, int d1) const
      {
        return 2*(d0+d1);
      }

      static constexpr result_type eval(T1 const& v1, T2 const& v2)
      {
        return std::atan2(v1, v2);
      }
      constexpr result_type operator()(T1 const& v1, T2 const& v2) const
      {
        return eval(v1, v2);
      }
    };

  } // end namespace functors
} // end namespace AMDiS
