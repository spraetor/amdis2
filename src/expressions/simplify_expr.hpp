/** \file simplify_expr.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
#include "value_expr.hpp"

#define SINGLE_ARG(...) __VA_ARGS__

namespace AMDiS 
{
  namespace traits 
  {
    template <class T>
    struct is_ct_value : false_ {};
    
    template <int I>
    struct is_ct_value< expressions::CValue<I> > : true_ {};
    
  } // end namespace traits
  
    
  namespace expressions 
  {
    template <class Term, class Enabled = void>
    struct Simplify 
    { 
      typedef Term type; 
      
      template <class Term_>
      static type eval(Term_&& t) { return t; }
    }; // by default do not simplify a term
    
  } // end namespace expressions
  

  // generator function for simplification
  template <class Term> 
  inline typename expressions::Simplify<Term>::type 
  simplify(Term&& t) 
  { return expressions::Simplify<Term>::eval(std::forward<Term>(t)); }


  namespace expressions 
  {
    /// -(N) -> (-N)
    template <int N>
    struct Simplify< Negative<CValue<N> > > 
    { 
      typedef CValue<-N> type; 
      
      static type eval(Negative<CValue<N> > const& t)
      { return {}; }
      
    };
    
    /// -0 -> 0
    template <>
    struct Simplify< Negative<CValue<0> > > 
    { 
      typedef CValue<0> type;
      
      static type eval(Negative<CValue<0> > const& t) 
      { return {}; }
    };
    
    /// (N) + (M) -> (N+M)
    template <int N, int M>
    struct Simplify< Add<CValue<N>, CValue<M> > > 
    { 
      typedef CValue<N+M> type;
      
      static type eval(Add<CValue<N>, CValue<M> > const& t) 
      { return {}; }
    };
    
    /// (N) * (M) -> (N*M)
    template <int N, int M >
    struct Simplify< Mult<CValue<N>, CValue<M> > > 
    { 
      typedef CValue<N*M> type;
      
      static type eval(Mult<CValue<N>, CValue<M> > const& t) 
      { return {}; }
    };
    
    /// (M) ^ N -> (M^N)
    template <int N, int M >
    struct Simplify< Pow<N, CValue<M> > > 
    { 
      typedef CValue<meta::pow<M,N>::value> type;
      
      static type eval(Pow<N, CValue<M> > const& t) 
      { return {}; }
    };
    
  } // end namespace expressions

  namespace expressions 
  {
    /// X + 0 -> X
    template <class Term>
    struct Simplify< Add<Term, CValue<0> >, typename enable_if_c<!traits::is_ct_value<Term>::value>::type >
    { 
      typedef typename Simplify<Term>::type type; 
      
      static type eval(Add<Term, CValue<0> > const& t) 
      { return simplify(t.getTerm(_0)); }
    };
    
    /// 0 + X -> X
    template <class Term>
    struct Simplify< Add<CValue<0>, Term>, typename enable_if_c<!traits::is_ct_value<Term>::value>::type >
    { 
      typedef typename Simplify<Term>::type type; 
      
      static type eval(Add<CValue<0>, Term> const& t) 
      { return simplify(t.getTerm(_1)); }
    };
      
    /// X * 0 -> 0
    template <class Term>
    struct Simplify< Mult<Term, CValue<0> >, typename enable_if_c<!traits::is_ct_value<Term>::value>::type > 
    { 
      typedef CValue<0> type; 
      
      static type eval(Mult<Term, CValue<0> > const& t) 
      { return {}; }
    };
    
    /// 0 * X -> 0
    template <class Term>
    struct Simplify< Mult<CValue<0>, Term>, typename enable_if_c<!traits::is_ct_value<Term>::value>::type > 
    { 
      typedef CValue<0> type; 
      
      static type eval(Mult<CValue<0>, Term> const& t) 
      { return {}; }
    };
      
    /// X * 1 -> X
    template <class Term>
    struct Simplify< Mult<Term, CValue<1> >, typename enable_if_c<!traits::is_ct_value<Term>::value>::type > 
    { 
      typedef typename Simplify<Term>::type type; 
      
      static type eval(Mult<Term, CValue<1> > const& t) 
      { return simplify(t.getTerm(_0)); }
    };
    
    /// 1 * X -> X
    template <class Term>
    struct Simplify< Mult<CValue<1>, Term>, typename enable_if_c<!traits::is_ct_value<Term>::value>::type > 
    { 
      typedef typename Simplify<Term>::type type; 
      
      static type eval(Mult<CValue<1>, Term> const& t) 
      { return simplify(t.getTerm(_1)); }
    };
    
  } // end namespace expressions

    
  namespace expressions 
  {
    /// X / X -> 1
    template <class Term>
    struct Simplify< Mult<Term, MultInverse<Term> > > 
    { 
      typedef CValue<1> type;
      
      static type eval(Mult<Term, MultInverse<Term> > const& t) 
      { return {}; }    
    };
    
    /// X / X -> 1
    template <class Term>
    struct Simplify< Mult<MultInverse<Term>, Term> > 
    { 
      typedef CValue<1> type;
      
      static type eval(Mult<MultInverse<Term>, Term> const& t)  
      { return {}; }
    };
    
    /// -(-X) -> X
    template <class Term>
    struct Simplify< Negative<Negative<Term> > > 
    { 
      typedef Term type;
      
      static type eval(Negative<Negative<Term> > const& t)  
      { return t.getTerm(_0).getTerm(_0); }
    };
    
    /// (X^(-1))^(-1) -> X
    template <class Term>
    struct Simplify< MultInverse<MultInverse<Term> > > 
    { 
      typedef Term type;
      
      static type eval(MultInverse<MultInverse<Term> > const& t)  
      { return t.getTerm(_0).getTerm(_0); }
    };
    
    /// -(X - Y) -> Y - X
    template <class X, class Y>
    struct Simplify< Negative<Add<X, Negative<Y> > > > 
    { 
      typedef Subtract<Y, X> type;
      
      static type eval(Negative<Add<X, Negative<Y> > > const& t) 
      { return t.getTerm(_0).getTerm(_1).getTerm(_0) - t.getTerm(_0).getTerm(_0); }
    };
    
    /// X + (-Y) -> X - Y
    template <class X, class Y>
    struct Simplify< Add<X, Negative<Y> > > 
    { 
      typedef Subtract<X, Y> type;
      
      static type eval(Add<X, Negative<Y> > const& t) 
      { return t.getTerm(_0) - t.getTerm(_1).getTerm(_0); }
    };
    
  } // end namespace expressions
  
    
  namespace expressions 
  {
    /// X^1 -> X
    template <class Term>
    struct Simplify< Pow<1, Term>, typename enable_if_c<!traits::is_ct_value<Term>::value>::type >
    { 
      typedef Term type; 
      
      static type eval(Pow<1,Term> const& t) { return t.getTerm(_0); }
    };
    
    /// X^0 -> 1
    template <class Term>
    struct Simplify< Pow<0, Term>, typename enable_if_c<!traits::is_ct_value<Term>::value>::type >
    { 
      typedef CValue<1> type;
      
      static type eval(Pow<0,Term> const& t) { return {}; }
    };

  } // end namespace expressions
  
    
  namespace expressions 
  {
    /// (XY) / (XZ) -> Y/Z
    template <class X, class Y, class Z>
    struct Simplify< Mult<Mult<X,Y>, MultInverse<Mult<X,Z> > > > 
    { 
      typedef Mult<Y, MultInverse<Z> > type; 
      
      // Y = t.getTerm(_0).getTerm(_1), Z = t.getTerm(_1).term.getTerm(_1)
      static type eval(Mult<Mult<X,Y>, MultInverse<Mult<X,Z> > > const& t) 
      { return t.getTerm(_0).getTerm(_1) / t.getTerm(_1).getTerm(_0).getTerm(_1); }    
    };
    
    /// (XY) / X -> Y
    template <class X, class Y>
    struct Simplify< Mult<Mult<X,Y>, MultInverse<X> > > 
    { 
      typedef Y type;
      
      // Y = t.getTerm(_0).getTerm(_1)
      static type eval(Mult<Mult<X,Y>, MultInverse<X> > const& t)  
      { return t.getTerm(_0).getTerm(_1); }
    };
    
    /// X / (XY) -> 1/Y
    template <class X, class Y>
    struct Simplify< Mult<X, MultInverse<Mult<X,Y> > > > 
    { 
      typedef MultInverse<Y> type;
      
      // Y = t.getTerm(_1).term.getTerm(_1)
      static type eval(Mult<Mult<X,Y>, MultInverse<X> > const& t)   
      { return MultInverse<Y>(t.getTerm(_1).getTerm(_0).getTerm(_1)); }
    };
    
    /// N*(M*X) -> (N*M) * X
    template <int N, int M, class X>
    struct Simplify< Mult<CValue<N>, Mult<CValue<M>, X> >, typename enable_if_c< N!=0 && N!=1 && M!=0 && M!=1 >::type  > 
    { 
      typedef Mult<CValue<N*M>, X> type;
      
      // X = t.getTerm(_1).getTerm(_1)
      static type eval(Mult<CValue<N>, Mult<CValue<M>, X> > const& t)    
      { return CValue<N*M>() * (t.getTerm(_1).getTerm(_1)); } 
    };
      
    /// N*(X*M) -> (N*M) * X
    template <int N, int M, class X>
    struct Simplify< Mult<CValue<N>, Mult<X, CValue<M> > >, typename enable_if_c< N!=0 && N!=1 && M!=0 && M!=1 >::type > 
    { 
      typedef Mult<CValue<N*M>, X> type; 
      
      // X = t.getTerm(_1).getTerm(_0)
      static type eval(Mult<CValue<N>, Mult<X, CValue<M> > > const& t)    
      { return CValue<N*M>() * (t.getTerm(_1).getTerm(_0)); } 
    };
    
    
    /// (M*X)*N -> (N*M) * X
    template <int N, int M, class X>
    struct Simplify< Mult<Mult<CValue<M>, X>, CValue<N> >, typename enable_if_c< N!=0 && N!=1 && M!=0 && M!=1 >::type > 
    { 
      typedef Mult<CValue<N*M>, X> type; 
      
      // X = t.getTerm(_0).getTerm(_1)
      static type eval(Mult<Mult<CValue<M>, X>, CValue<N> > const& t)    
      { return CValue<N*M>() * (t.getTerm(_0).getTerm(_1)); }
    };
    
    
    /// (X*M)*N -> (N*M) * X
    template <int N, int M, class X>
    struct Simplify< Mult<Mult<X, CValue<M> >, CValue<N> >, typename enable_if_c< N!=0 && N!=1 && M!=0 && M!=1 >::type > 
    { 
      typedef Mult<CValue<N*M>, X> type; 
      
      // X = t.getTerm(_0).getTerm(_0)
      static type eval(Mult<Mult<X, CValue<M> >, CValue<N> > const& t)    
      { return CValue<N*M>() * (t.getTerm(_0).getTerm(_0)); }
    };

  } // end namespace expressions


  namespace expressions 
  {
    /// A + B -> simplify(A) + simplify(B)
    template <class Term1, class Term2>
    class Simplify< Add<Term1, Term2>, typename enable_if_c<!(traits::is_ct_value<Term1>::value || traits::is_ct_value<Term2>::value)>::type > 
    { 
      typedef typename Simplify<Term1>::type S1;
      typedef typename Simplify<Term2>::type S2;
      
    public:
      typedef Add<S1, S2> type;
      
      static type eval(Add<Term1, Term2> const& t) 
      {
        return simplify(t.getTerm(_0)) + simplify(t.getTerm(_1));
      }
    };
    
    /// A * B -> simplify(A) * simplify(B)
    template <class Term1, class Term2>
    class Simplify< Mult<Term1, Term2>, typename enable_if_c<!(traits::is_ct_value<Term1>::value || traits::is_ct_value<Term2>::value)>::type > 
    { 
      typedef typename Simplify<Term1>::type S1;
      typedef typename Simplify<Term2>::type S2;
      
    public:
      typedef Mult<S1, S2> type;
      
      static type eval(Mult<Term1, Term2> const& t) 
      {
        return simplify(t.getTerm(_0)) * simplify(t.getTerm(_1));
      }
    };
    
    /// X(A) -> X(simplify(A))
    template <template<class> class Outer, class Inner>
    class Simplify< Outer<Inner>, typename enable_if_c<!(traits::is_ct_value<Inner>::value) && traits::is_expr<Inner>::value>::type > 
    { 
      typedef typename Simplify<Inner>::type S;
      
    public:
      typedef Outer<S> type;
      
      static type eval(Outer<Inner> const& t) 
      {
        return type(simplify(t.getTerm(_0)));
      }
    };

  } // end namespace expressions


  // =============================================================================
  // multiple simplifications

  namespace expressions 
  {    
    template <int N, class Term, class Enabled = void>
    struct SimplifyRecursive 
    { 
      typedef typename Simplify< typename SimplifyRecursive<N-1, Term>::type >::type type;
    };
    
    template <class Term>
    struct SimplifyRecursive<1, Term>
    { 
      typedef typename Simplify< Term >::type type;
    };
    
    template <int N, class Term >
    struct SimplifyRecursive<N, Term, typename enable_if_c<(N <= 0)>::type>
    { 
      typedef Term type;
    };
    
  } // end namespace expressions

    
  // forward declaration
  template <int N, class Term> 
  inline typename expressions::SimplifyRecursive<N, Term>::type
  simplify(const Term& t);


  template <int N, class Term> 
  inline typename expressions::SimplifyRecursive<N, Term>::type
  simplify(Term&& t, int_<N>) { return simplify(simplify<N-1>(std::forward<Term>(t))); }

  template <class Term> 
  inline typename expressions::SimplifyRecursive<1, Term>::type
  simplify(Term&& t, int_<1>) { return simplify(std::forward<Term>(t)); }

  template <class Term> 
  inline typename expressions::SimplifyRecursive<0, Term>::type
  simplify(Term&& t, int_<0>) { return t; }


  template <int N, class Term> 
  inline typename expressions::SimplifyRecursive<N, Term>::type
  simplify(Term&& t) { return simplify(std::forward<Term>(t), int_<N>()); }

} // end namespace AMDiS

#endif // AMDIS_SIMPLIFY_EXPRESSION_HPP
