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



/** \file meta.hpp */

#pragma once

#include <boost/mpl/logical.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/not.hpp>

#include <boost/version.hpp> 
#if BOOST_VERSION >= 105600
  #include <boost/core/enable_if.hpp>
#else
  #include <boost/utility/enable_if.hpp>
#endif

namespace AMDiS 
{
  // introduce some shortcuts for boost::mpl
  // ---------------------------------------
  using boost::mpl::bool_;
  using boost::mpl::true_;
  using boost::mpl::false_;
  using boost::mpl::and_;
  using boost::mpl::or_;
  using boost::mpl::not_;
  using boost::mpl::int_;
  
  using boost::mpl::if_;
  using boost::mpl::if_c;
  
  using boost::enable_if;
  using boost::disable_if;
  
  using boost::enable_if_c;
  using boost::disable_if_c;
  
  /// Namespace for compile-time functions
  /** The namespace \ref meta provides some meta-functions mainly for
    *  compile-time integer arithmetics
    **/
  namespace meta { }
  
  // _________________________________________________________________________ 
#if HAS_CONSTEXPR
    /// MAX(x,y)  
    template <class T> constexpr T MAX(T x, T y) { return x > y ? x : y; }
#else
    namespace meta {
      /// max<X, Y>   
      template <int X, int Y> struct max : int_< (X > Y ? X : Y) > {};
    }
    #define MAX(X,Y) meta::max<X,Y>::value
#endif
    
  // _________________________________________________________________________
#if HAS_CONSTEXPR
    /// MIN(x,y)
    template <class T> constexpr T MIN(T x, T y) { return x < y ? x : y; }
#else
    namespace meta {
      /// min<X, Y>
      template <int X, int Y> struct min : int_< (X > Y ? Y : X) > {};
    }
    #define MIN(X,Y) min<X,Y>::value
#endif
    
  namespace meta
  {
    // _________________________________________________________________________
    /// abs<X> == |X|
    template <int X>
    struct abs : int_< (X >= 0 ? X : -X) > {};

    // _________________________________________________________________________
    /// pow<X, p> == X^p
    template <int X, int P>
    struct pow : int_< (X * pow<X, P-1>::value) > {};

    /// \cond HIDDEN_SYMBOLS
    template <int X>
    struct pow<X, 0> : int_< 1 > {};
    /// \endcond

    /// sqr<X> = X^2
    template <int X>
    struct sqr : pow<X, 2> {};
  
    // _________________________________________________________________________
    /// root<x, p> == p-th-root(x)
    template <int X, int P, int I = 1>
    struct root : if_c
      < 
        boost::mpl::less<pow<I, P>, int_<X> >::value, 
        root<X, P, I+1>, 
        int_< I > 
      >::type {};

    /// \cond HIDDEN_SYMBOLS
    template <int X, int P>
    struct root<X, P, X> : int_< X > {};
    /// \endcond

    /// sqrt<X> = integer sqrt of X
    template <int X>
    struct sqrt : root<X, 2> {};

    // _________________________________________________________________________
    /// log<x, p> == log_p(x)
    template <int X, int P, int I = 1>
    struct log : if_c
      < (pow<P, I>::value < X), log<X, P, I+1>, int_< I > >::type {};

    /// \cond HIDDEN_SYMBOLS
    template <int X, int P>
    struct log<X, P, X> : int_< 0 > {};
    /// \endcond

    /// evaluates to true if \p X is power of \p Base
    template <int X, int Base>
    struct is_power_of : bool_< (pow<Base, log<X, Base>::value>::value == X) > {};

    // _________________________________________________________________________
    /// evalutes to true, if \p U is divisible by \p V
    template <int U, int V>
    struct is_divisible_by : bool_< (U % V == 0) > {};

    /// \cond HIDDEN_SYMBOLS
    template <int U>
    struct is_divisible_by<U, 0> : false_ {};
    /// \endcond

    /// evalutes to true, if \p U is divisible by 2
    template <int U>
    struct is_even : is_divisible_by<U, 2> {};
    
    /// evalutes to true, if \p U is not divisible by 2
    template <int U>
    struct is_odd : not_< is_even<U> > {};

    
  } // end namespace meta

} // end namespace AMDiS
