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

#ifndef AMDIS_OPERATIONS_META_HPP
#define AMDIS_OPERATIONS_META_HPP

#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/not.hpp>

namespace AMDiS 
{
  namespace meta
  {
    /// abs<X> == |X|
    // _________________________________________________________________________
    template<int X>
    struct abs : boost::mpl::int_< (X >= 0 ? X : -X) > {};

    /// pow<x, p> == x^p
    // _________________________________________________________________________
    template<int X, int P>
    struct pow : boost::mpl::int_< (X * pow<X, P-1>::value) > {};

    template<int X>
    struct pow<X, 0> : boost::mpl::int_< 1 > {};

    template<int X>
    struct sqr : pow<X, 2> {};
  
    /// root<x, p> == p-th-root(x)
    // _________________________________________________________________________
    template<int X, int P, int I = 1>
    struct root : boost::mpl::if_c
      < 
        boost::mpl::less<pow<I, P>, boost::mpl::int_<X> >::value, 
        root<X, P, I+1>, 
        boost::mpl::int_< I > 
      >::type {};

    template<int X, int P>
    struct root<X, P, X> : boost::mpl::int_< X > {};

    template<int X>
    struct sqrt : root<X, 2> {};

    /// log<x, p> == log_p(x)
    // _________________________________________________________________________
    template<int X, int P, int I = 1>
    struct log : boost::mpl::if_c
      < 
        boost::mpl::less<pow<P, I>, boost::mpl::int_<X> >::value, 
        log<X, P, I+1>, 
        boost::mpl::int_< I > 
      >::type {};

    template<int X, int P>
    struct log<X, P, X> : boost::mpl::int_< 0 > {};


    template<int X, int Base>
    struct is_power_of : boost::mpl::bool_< (pow<Base, log<X, Base>::value>::value == X) > {};

    /// is divisible by
    // _________________________________________________________________________
    template<int U, int V>
    struct is_divisible_by : boost::mpl::bool_< (U % V == 0) > {};

    template<int U>
    struct is_divisible_by<U, 0> : boost::mpl::false_ {};

    template<int U>
    struct is_even : is_divisible_by<U, 2> {};

    template<int U>
    struct is_odd : boost::mpl::not_< is_even<U> > {};

    /// generic loops
    // _________________________________________________________________________
    template<long I, long N> 
    struct FOR 
    { 
      template<class A, class B, class Op>
      static void transform(A const& a, B& b, Op op, size_t shift = 0)
      {
	b[I+shift] = op(a[I+shift]);
	FOR<I+1,N>::transform(a, b, op, shift);
      }
      
      template<class A, class T, class Op, class BinaryOp>
      static T accumulate(A const& a, T init, Op op, BinaryOp binary_op, size_t shift = 0)
      {
	return binary_op(FOR<I+1,N>::accumulate(a, init, op, binary_op, shift), 
			 op(a[I+shift]));
      }
      
      template<class A, class B, class T, class BinaryOp1, class BinaryOp2>
      static T inner_product(A const& a, B const& b, T init, BinaryOp1 binary_op1, BinaryOp2 binary_op2, size_t shift = 0)
      {
	return binary_op1(FOR<I+1,N>::inner_product(a, b, init, binary_op1, binary_op2, shift), 
			  binary_op2(a[I+shift], b[I+shift]));
      }
    };

    // Abbruchbedingung I==N
    template<long N> 
    struct FOR<N, N>
    { 
      template<class A, class B, class Op>
      static void transform(A const& a, B& b, Op op, size_t shift = 0) {}
      
      template<class A, class T, class Op, class BinaryOp>
      static T accumulate(A const& a, T init, Op op, BinaryOp binary_op, size_t shift = 0)
      {
	return init;
      }
      
      template<class A, class B, class T, class BinaryOp1, class BinaryOp2>
      static T inner_product(A const& a, B const& b, T init, BinaryOp1 binary_op1, BinaryOp2 binary_op2, size_t shift = 0)
      {
	return init;
      }
    };
    
  } // end namespace meta

} // end namespace AMDiS

#endif // AMDIS_OPERATIONS_META_HPP
