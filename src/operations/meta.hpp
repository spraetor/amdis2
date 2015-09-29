/** \file meta.hpp */

#pragma once

#include <traits/meta_basic.hpp>

namespace AMDiS
{
  namespace meta
  {
    /// abs<X> == |X|
    // _________________________________________________________________________
    template <int X>
    using abs = int_< (X >= 0 ? X : -X) >;

    /// pow<x, p> == x^p
    // _________________________________________________________________________
    template <int X, int P>
    struct pow : int_< (X * pow<X, P-1>::value) > {};

    template <int X>
    struct pow<X, 0> : int_<1> {};

    template <int X>
    using sqr = pow<X, 2>;

    /// root<x, p> == p-th-root(x)
    // _________________________________________________________________________
    template <int X, int P, int I = 1>
    struct root : if_then_else
      <
    (pow<I, P>::value < X),
    root<X, P, I+1>,
         int_<I>
         >::type {};

    template <int X, int P>
    struct root<X, P, X> : int_<X> {};

    template <int X>
    using sqrt = root<X, 2>;

    /// log<x, p> == log_p(x)
    // _________________________________________________________________________
    template <int X, int P, int I = 1>
    struct log : if_then_else
      <
    (pow<P, I>::value < X),
    log<X, P, I+1>,
        int_<I>
        >::type {};

    template <int X, int P>
    struct log<X, P, X> : int_<0> {};


    template <int X, int Base>
    using is_power_of = bool_<(pow<Base, log<X, Base>::value>::value == X)>;

    /// is divisible by
    // _________________________________________________________________________
    template <int U, int V>
    struct is_divisible_by : bool_< (U % V == 0) > {};

    template <int U>
    struct is_divisible_by<U, 0> : false_ {};

    template <int U>
    using is_even = is_divisible_by<U, 2>;

    template <int U>
    using is_odd = not_<is_even<U>>;

  } // end namespace meta

} // end namespace AMDiS
