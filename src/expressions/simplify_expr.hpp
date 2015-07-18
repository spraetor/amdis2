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



/** \file simplify_expr.hpp */

#ifndef AMDIS_SIMPLIFY_EXPRESSION_HPP
#define AMDIS_SIMPLIFY_EXPRESSION_HPP

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
#include "value_expr.hpp"

#define SINGLE_ARG(...) __VA_ARGS__

namespace AMDiS 
{
  namespace traits 
  {
    template<typename T>
    struct is_ct_value : boost::mpl::false_ {};
    
    template<int I>
    struct is_ct_value< expressions::CValue<I> > : boost::mpl::true_ {};
    
  } // end namespace traits
  
    
  namespace expressions 
  {
    template<typename Term, typename Enabled = void>
    struct Simplify 
    { 
      typedef Term type; 
      
      static type eval(Term const& t) { return t; }
    }; // by default do not simplify a term
    
  } // end namespace expressions
  

  // generator function for simplification
  template<typename Term> 
  inline typename expressions::Simplify<Term>::type simplify(const Term& t) 
  { return expressions::Simplify<Term>::eval(t); }


  namespace expressions 
  {
    /// -(N) -> (-N)
    template<int N>
    struct Simplify< Negative<CValue<N> > > 
    { 
      typedef CValue<-N> type; 
      
      static type eval(Negative<CValue<N> > const& t)
      { return CValue<-N>(); }
      
    };
    
    /// -0 -> 0
    template<>
    struct Simplify< Negative<CValue<0> > > 
    { 
      typedef CValue<0> type;
      
      static type eval(Negative<CValue<0> > const& t) 
      { return CValue<0>(); }
    };
    
    /// (N) + (M) -> (N+M)
    template<int N, int M>
    struct Simplify< Add<CValue<N>, CValue<M> > > 
    { 
      typedef CValue<N+M> type;
      
      static type eval(Add<CValue<N>, CValue<M> > const& t) 
      { return CValue<N+M>(); }
    };
    
    /// (N) * (M) -> (N*M)
    template<int N, int M >
    struct Simplify< Mult<CValue<N>, CValue<M> > > 
    { 
      typedef CValue<N*M> type;
      
      static type eval(Mult<CValue<N>, CValue<M> > const& t) 
      { return CValue<N*M>(); }
    };
    
    /// (M) ^ N -> (M^N)
    template<int N, int M >
    struct Simplify< Pow<N, CValue<M> > > 
    { 
      typedef CValue<meta::pow<M,N>::value> type;
      
      static type eval(Pow<N, CValue<M> > const& t) 
      { return type(); }
    };
    
  } // end namespace expressions

  namespace expressions 
  {
    /// X + 0 -> X
    template<typename Term>
    struct Simplify< Add<Term, CValue<0> >, typename boost::enable_if_c<!traits::is_ct_value<Term>::value>::type >
    { 
      typedef typename Simplify<Term>::type type; 
      
      static type eval(Add<Term, CValue<0> > const& t) 
      { return simplify(t.term1); }
    };
    
    /// 0 + X -> X
    template<typename Term>
    struct Simplify< Add<CValue<0>, Term>, typename boost::enable_if_c<!traits::is_ct_value<Term>::value>::type >
    { 
      typedef typename Simplify<Term>::type type; 
      
      static type eval(Add<CValue<0>, Term> const& t) 
      { return simplify(t.term2); }
    };
      
    /// X * 0 -> 0
    template<typename Term>
    struct Simplify< Mult<Term, CValue<0> >, typename boost::enable_if_c<!traits::is_ct_value<Term>::value>::type > 
    { 
      typedef CValue<0> type; 
      
      static type eval(Mult<Term, CValue<0> > const& t) 
      { return CValue<0>(); }
    };
    
    /// 0 * X -> 0
    template<typename Term>
    struct Simplify< Mult<CValue<0>, Term>, typename boost::enable_if_c<!traits::is_ct_value<Term>::value>::type > 
    { 
      typedef CValue<0> type; 
      
      static type eval(Mult<CValue<0>, Term> const& t) 
      { return CValue<0>(); }
    };
      
    /// X * 1 -> X
    template<typename Term>
    struct Simplify< Mult<Term, CValue<1> >, typename boost::enable_if_c<!traits::is_ct_value<Term>::value>::type > 
    { 
      typedef typename Simplify<Term>::type type; 
      
      static type eval(Mult<Term, CValue<1> > const& t) 
      { return simplify(t.term1); }
    };
    
    /// 1 * X -> X
    template<typename Term>
    struct Simplify< Mult<CValue<1>, Term>, typename boost::enable_if_c<!traits::is_ct_value<Term>::value>::type > 
    { 
      typedef typename Simplify<Term>::type type; 
      
      static type eval(Mult<CValue<1>, Term> const& t) 
      { return simplify(t.term2); }
    };
    
  } // end namespace expressions

    
  namespace expressions 
  {
    /// X / X -> 1
    template<typename Term>
    struct Simplify< Mult<Term, MultInverse<Term> > > 
    { 
      typedef CValue<1> type;
      
      static type eval(Mult<Term, MultInverse<Term> > const& t) 
      { return CValue<1>(); }    
    };
    
    /// X / X -> 1
    template<typename Term>
    struct Simplify< Mult<MultInverse<Term>, Term> > 
    { 
      typedef CValue<1> type;
      
      static type eval(Mult<MultInverse<Term>, Term> const& t)  
      { return CValue<1>(); }
    };
    
    /// -(-X) -> X
    template<typename Term>
    struct Simplify< Negative<Negative<Term> > > 
    { 
      typedef Term type;
      
      static type eval(Negative<Negative<Term> > const& t)  
      { return t.term.term; }
    };
    
    /// (X^(-1))^(-1) -> X
    template<typename Term>
    struct Simplify< MultInverse<MultInverse<Term> > > 
    { 
      typedef Term type;
      
      static type eval(MultInverse<MultInverse<Term> > const& t)  
      { return t.term.term; }
    };
    
    /// -(X - Y) -> Y - X
    template<typename X, typename Y>
    struct Simplify< Negative<Add<X, Negative<Y> > > > 
    { 
      typedef Add<Y, Negative<X> > type;
      
      static type eval(Negative<Add<X, Negative<Y> > > const& t) 
      { return t.term.term2.term - t.term.term1; }
    };
    
  } // end namespace expressions
  
    
  namespace expressions 
  {
    /// X^1 -> X
    template<typename Term>
    struct Simplify< Pow<1, Term>, typename boost::enable_if_c<!traits::is_ct_value<Term>::value>::type >
    { 
      typedef Term type; 
      
      static type eval(Pow<1,Term> const& t) { return t.term; }
    };
    
    /// X^0 -> 1
    template<typename Term>
    struct Simplify< Pow<0, Term>, typename boost::enable_if_c<!traits::is_ct_value<Term>::value>::type >
    { 
      typedef CValue<1> type;
      
      static type eval(Pow<0,Term> const& t) { return CValue<1>(); }
    };

  } // end namespace expressions
  
    
  namespace expressions 
  {
    /// (XY) / (XZ) -> Y/Z
    template<typename X, typename Y, typename Z>
    struct Simplify< Mult<Mult<X,Y>, MultInverse<Mult<X,Z> > > > 
    { 
      typedef Mult<Y, MultInverse<Z> > type; 
      
      // Y = t.term1.term2, Z = t.term2.term.term2
      static type eval(Mult<Mult<X,Y>, MultInverse<Mult<X,Z> > > const& t) 
      { return t.term1.term2 / t.term2.term.term2; }    
    };
    
    /// (XY) / X -> Y
    template<typename X, typename Y>
    struct Simplify< Mult<Mult<X,Y>, MultInverse<X> > > 
    { 
      typedef Y type;
      
      // Y = t.term1.term2
      static type eval(Mult<Mult<X,Y>, MultInverse<X> > const& t)  
      { return t.term1.term2; }
    };
    
    /// X / (XY) -> 1/Y
    template<typename X, typename Y>
    struct Simplify< Mult<X, MultInverse<Mult<X,Y> > > > 
    { 
      typedef MultInverse<Y> type;
      
      // Y = t.term2.term.term2
      static type eval(Mult<Mult<X,Y>, MultInverse<X> > const& t)   
      { return MultInverse<Y>(t.term2.term.term2); }
    };
    
    /// N*(M*X) -> (N*M) * X
    template<int N, int M, typename X>
    struct Simplify< Mult<CValue<N>, Mult<CValue<M>, X> >, typename boost::enable_if_c< N!=0 && N!=1 && M!=0 && M!=1 >::type  > 
    { 
      typedef Mult<CValue<N*M>, X> type;
      
      // X = t.term2.term2
      static type eval(Mult<CValue<N>, Mult<CValue<M>, X> > const& t)    
      { return CValue<N*M>() * (t.term2.term2); } 
    };
      
    /// N*(X*M) -> (N*M) * X
    template<int N, int M, typename X>
    struct Simplify< Mult<CValue<N>, Mult<X, CValue<M> > >, typename boost::enable_if_c< N!=0 && N!=1 && M!=0 && M!=1 >::type > 
    { 
      typedef Mult<CValue<N*M>, X> type; 
      
      // X = t.term2.term1
      static type eval(Mult<CValue<N>, Mult<X, CValue<M> > > const& t)    
      { return CValue<N*M>() * (t.term2.term1); } 
    };
    
    
    /// (M*X)*N -> (N*M) * X
    template<int N, int M, typename X>
    struct Simplify< Mult<Mult<CValue<M>, X>, CValue<N> >, typename boost::enable_if_c< N!=0 && N!=1 && M!=0 && M!=1 >::type > 
    { 
      typedef Mult<CValue<N*M>, X> type; 
      
      // X = t.term1.term2
      static type eval(Mult<Mult<CValue<M>, X>, CValue<N> > const& t)    
      { return CValue<N*M>() * (t.term1.term2); }
    };
    
    
    /// (X*M)*N -> (N*M) * X
    template<int N, int M, typename X>
    struct Simplify< Mult<Mult<X, CValue<M> >, CValue<N> >, typename boost::enable_if_c< N!=0 && N!=1 && M!=0 && M!=1 >::type > 
    { 
      typedef Mult<CValue<N*M>, X> type; 
      
      // X = t.term1.term1
      static type eval(Mult<Mult<X, CValue<M> >, CValue<N> > const& t)    
      { return CValue<N*M>() * (t.term1.term1); }
    };

  } // end namespace expressions


  namespace expressions 
  {
    /// A + B -> simplify(A) + simplify(B)
    template<typename Term1, typename Term2>
    class Simplify< Add<Term1, Term2>, typename boost::enable_if_c<!(traits::is_ct_value<Term1>::value || traits::is_ct_value<Term2>::value)>::type > 
    { 
      typedef typename Simplify<Term1>::type S1;
      typedef typename Simplify<Term2>::type S2;
      
    public:
      typedef Add<S1, S2> type;
      
      static type eval(Add<Term1, Term2> const& t) 
      {
	return simplify(t.term1) + simplify(t.term2);
      }
    };
    
    /// A * B -> simplify(A) * simplify(B)
    template<typename Term1, typename Term2>
    class Simplify< Mult<Term1, Term2>, typename boost::enable_if_c<!(traits::is_ct_value<Term1>::value || traits::is_ct_value<Term2>::value)>::type > 
    { 
      typedef typename Simplify<Term1>::type S1;
      typedef typename Simplify<Term2>::type S2;
      
    public:
      typedef Mult<S1, S2> type;
      
      static type eval(Mult<Term1, Term2> const& t) 
      {
	return simplify(t.term1) * simplify(t.term2);
      }
    };
    
    /// X(A) -> X(simplify(A))
    template<template<class> class Outer, class Inner>
    class Simplify< Outer<Inner>, typename boost::enable_if_c<!(traits::is_ct_value<Inner>::value) && traits::is_expr<Inner>::value>::type > 
    { 
      typedef typename Simplify<Inner>::type S;
      
    public:
      typedef Outer<S> type;
      
      static type eval(Outer<Inner> const& t) 
      {
	return type(simplify(t.term));
      }
    };

  } // end namespace expressions


  // =============================================================================
  // multiple simplifications

  namespace expressions 
  {    
    template<int N, typename Term, typename Enabled = void>
    struct SimplifyRecursive 
    { 
      typedef typename Simplify< typename SimplifyRecursive<N-1, Term>::type >::type type;
    };
    
    template<typename Term>
    struct SimplifyRecursive<1, Term>
    { 
      typedef typename Simplify< Term >::type type;
    };
    
    template<int N, typename Term >
    struct SimplifyRecursive<N, Term, typename boost::enable_if_c<(N <= 0)>::type>
    { 
      typedef Term type;
    };
    
  } // end namespace expressions

    
  template<int N, typename Term> 
  inline typename expressions::SimplifyRecursive<N, Term>::type
  simplify(const Term& t);

  template<int N, typename Term> 
  inline typename expressions::SimplifyRecursive<N, Term>::type
  simplify(const Term& t, boost::mpl::int_<N>) { return simplify(simplify<N-1>(t)); }

  template<typename Term> 
  inline typename expressions::SimplifyRecursive<1, Term>::type
  simplify(const Term& t, boost::mpl::int_<1>) { return simplify(t); }

  template<typename Term> 
  inline typename expressions::SimplifyRecursive<0, Term>::type
  simplify(const Term& t, boost::mpl::int_<0>) { return t; }

  template<int N, typename Term> 
  inline typename expressions::SimplifyRecursive<N, Term>::type
  simplify(const Term& t) { return simplify(t, boost::mpl::int_<N>()); }

} // end namespace AMDiS

#endif // AMDIS_SIMPLIFY_EXPRESSION_HPP
