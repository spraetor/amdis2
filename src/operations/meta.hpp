/** \file meta.hpp */

#pragma once

#include "traits/meta_basic.hpp"

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
    struct pow<X, 0> : int_< 1 > {};

    template <int X>
    using sqr = pow<X, 2>;
  
    /// root<x, p> == p-th-root(x)
    // _________________________________________________________________________
    template <int X, int P, int I = 1>
    struct root : if_then_else
      < 
        (pow<I, P>::value < X), 
        root<X, P, I+1>, 
        int_< I > 
      >::type {};

    template <int X, int P>
    struct root<X, P, X> : int_< X > {};

    template <int X>
    using sqrt = root<X, 2>;

    /// log<x, p> == log_p(x)
    // _________________________________________________________________________
    template <int X, int P, int I = 1>
    struct log : if_then_else
      < 
        (pow<P, I>::value < X), 
        log<X, P, I+1>, 
        int_< I > 
      >::type {};

    template <int X, int P>
    struct log<X, P, X> : int_< 0 > {};


    template <int X, int Base>
    using is_power_of = bool_< (pow<Base, log<X, Base>::value>::value == X) >;

    /// is divisible by
    // _________________________________________________________________________
    template <int U, int V>
    struct is_divisible_by : bool_< (U % V == 0) > {};

    template <int U>
    struct is_divisible_by<U, 0> : false_ {};

    template <int U>
    using is_even = is_divisible_by<U, 2>;

    template <int U>
    using is_odd = not_< is_even<U> >;

    /// generic loops
    // _________________________________________________________________________
    template <long I, long N> 
    struct FOR 
    { 
      template <class A, class B, class Op>
      static void transform(A const& a, B& b, Op op, size_t shift = 0)
      {
        b[I+shift] = op(a[I+shift]);
        FOR<I+1,N>::transform(a, b, op, shift);
      }
      
      template <class A, class T, class Op, class BinaryOp>
      static T accumulate(A const& a, T init, Op op, BinaryOp binary_op, size_t shift = 0)
      {
        return binary_op(FOR<I+1,N>::accumulate(a, init, op, binary_op, shift), 
                         op(a[I+shift]));
      }
      
      template <class A, class B, class T, class BinaryOp1, class BinaryOp2>
      static T inner_product(A const& a, B const& b, T init, BinaryOp1 binary_op1, BinaryOp2 binary_op2, size_t shift = 0)
      {
        return binary_op1(FOR<I+1,N>::inner_product(a, b, init, binary_op1, binary_op2, shift), 
                          binary_op2(a[I+shift], b[I+shift]));
      }
    };

    // Abbruchbedingung I==N
    template <long N> 
    struct FOR<N, N>
    { 
      template <class A, class B, class Op>
      static void transform(A const& a, B& b, Op op, size_t shift = 0) {}
      
      template <class A, class T, class Op, class BinaryOp>
      static T accumulate(A const& a, T init, Op op, BinaryOp binary_op, size_t shift = 0)
      {
        return init;
      }
      
      template <class A, class B, class T, class BinaryOp1, class BinaryOp2>
      static T inner_product(A const& a, B const& b, T init, BinaryOp1 binary_op1, BinaryOp2 binary_op2, size_t shift = 0)
      {
        return init;
      }
    };
    
  } // end namespace meta

} // end namespace AMDiS
