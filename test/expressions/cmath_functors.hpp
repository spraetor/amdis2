#pragma once

// std c++ headers
#include <cmath>
#include <type_traits>

// AMDiS includes
#include "Math.hpp"
#include "operations/functors.hpp"

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

    AMDIS_MAKE_UNARY_FUNCTOR( Ceil,    d0, std::ceil(v)  )
    AMDIS_MAKE_UNARY_FUNCTOR( Floor,   d0, std::floor(v) )
    AMDIS_MAKE_UNARY_FUNCTOR( Exp,   2*d0, std::exp(v)   )
    AMDIS_MAKE_UNARY_FUNCTOR( Log,   2*d0, std::log(v)   )
    AMDIS_MAKE_UNARY_FUNCTOR( Sin,   2*d0, std::sin(v)   )
    AMDIS_MAKE_UNARY_FUNCTOR( Cos,   2*d0, std::cos(v)   )
    AMDIS_MAKE_UNARY_FUNCTOR( Tan,   2*d0, std::tan(v)   )
    AMDIS_MAKE_UNARY_FUNCTOR( Asin,  2*d0, std::asin(v)  )
    AMDIS_MAKE_UNARY_FUNCTOR( Acos,  2*d0, std::acos(v)  )
    AMDIS_MAKE_UNARY_FUNCTOR( Atan,  2*d0, std::atan(v)  )
    AMDIS_MAKE_UNARY_FUNCTOR( Sinh,  2*d0, std::sinh(v)  )
    AMDIS_MAKE_UNARY_FUNCTOR( Cosh,  2*d0, std::cosh(v)  )
    AMDIS_MAKE_UNARY_FUNCTOR( Tanh,  2*d0, std::tanh(v)  )
    AMDIS_MAKE_UNARY_FUNCTOR( Asinh, 2*d0, std::asinh(v) )
    AMDIS_MAKE_UNARY_FUNCTOR( Acosh, 2*d0, std::acosh(v) )
    AMDIS_MAKE_UNARY_FUNCTOR( Atanh, 2*d0, std::atanh(v) )


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
