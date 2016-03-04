#pragma once

#include "traits/basic.hpp"
#include "traits/meta_basic.hpp"
#include "traits/scalar_types.hpp"
      
      
namespace AMDiS
{
  struct FunctorBase
  {
    int getDegree() const
    {
      return 0;
    }

    template <class... Ints>
      Requires_t<and_<concepts::Integral<Ints>...>, int>
    getDegree(Ints...) const
    {
      return 0;
    }
  };
} // end namespace AMDiS


/// Macro that generates a unary functor.
/**
 *  \p NAME    Name of the class.
 *  \p DEGREE  Expression in 'd0' that gives the polynomial degree, with
 *             'd0' the degree of the argument passed to the functor.
 *  \p FCT     Name of a unary c++-function that represents the functor.
 */
#ifndef AMDIS_MAKE_UNARY_FUNCTOR
#define AMDIS_MAKE_UNARY_FUNCTOR( NAME, DEGREE, FCT )     \
    template <class T>                                    \
    struct NAME : FunctorBase                             \
    {                                                     \
      using result_type = T;                              \
      constexpr int getDegree(int d0) const               \
      {                                                   \
        return DEGREE ;                                   \
      }                                                   \
      static constexpr auto eval(T const& v)       RETURNS( FCT )     \
      constexpr auto operator() (T const& v) const RETURNS( eval(v) ) \
    };
#endif

/// Macro that generates a binary functor.
/**
 *  \p NAME    Name of the class.
 *  \p DEGREE  Expression in 'd0,d1' that gives the polynomial degree, with
 *             'd0' the degree of the first argument and 'd1' the degree of
 *             the second argument passed to the functor.
 *  \p FCT     Expression in 'v0,v1' the represents the functor, with
 *             'v0' the first argument and 'v1' the second argument.
 */
#ifndef AMDIS_MAKE_BINARY_FUNCTOR
#define AMDIS_MAKE_BINARY_FUNCTOR( NAME, DEGREE, FCT )    \
    template <class T0, class T1 = T0>                    \
    struct NAME : FunctorBase                             \
    {                                                     \
      constexpr int getDegree(int d0, int d1) const       \
      {                                                   \
        return DEGREE ;                                   \
      }                                                   \
      static constexpr auto eval(T0 const& v0, T1 const& v1)       RETURNS( FCT )          \
      constexpr auto operator() (T0 const& v0, T1 const& v1) const RETURNS( eval(v0, v1) ) \
    };
#endif
