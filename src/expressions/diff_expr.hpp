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



/** \file diff_expr.hpp */

#ifndef AMDIS_DIFF_EXPRESSION_HPP
#define AMDIS_DIFF_EXPRESSION_HPP

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
#include "simplify_expr.hpp"

namespace AMDiS 
{    
  // by default the derivative is zero
  namespace expressions 
  {    
    struct DummyType {
      typedef _unknown id;
      typedef double value_type;
    };
    
    template<typename Identifier, typename Term, typename Direction = DummyType>
    struct Diff
    {
      typedef CValue<0> type;
      typedef CValue<0> dir_type;
      static type eval(Term const& t) { return CValue<0>(); }
      static dir_type eval(Term const& t, Direction const& d) { return CValue<0>(); }
    };
    
  } // end namespace expressions
    
  template<typename Id, typename Term>
  inline typename expressions::Diff<Id,Term>::type 
  diff(const Term& t) { return expressions::Diff<Id,Term>::eval(t); }
    
  template<typename Id, typename Term, typename Direction>
  inline typename expressions::Diff<Id,Term,Direction>::dir_type 
  diff(const Term& t, const Direction& d) { return expressions::Diff<Id,Term,Direction>::eval(t,d); }
  
  
  // _____________________________________________________________________________
  // value/gradient expressions
  namespace expressions 
  {      
    // ValueOf
    template<typename Id, typename Vector, typename Name, typename Direction>
    struct Diff< Id, ValueOf<Vector, Name>, Direction >
    {
      typedef ValueOf<Vector, Name> Term;
      typedef CValue<0> type;
      typedef CValue<0> dir_type;
      
      static type eval(Term const& t) { return CValue<0>(); }    
      static dir_type eval(Term const& t, Direction const& d) { return CValue<0>(); }
    };
    
    template<typename Id, typename Vector, typename Direction>
    struct Diff< Id, ValueOf<Vector, Id>, Direction >
    {
      typedef expressions::ValueOf<Vector, Id> Term;
      typedef CValue<1> type;
      typedef Direction dir_type;
      
      static type eval(Term const& t) { return CValue<1>(); }    
      static dir_type eval(Term const& t, Direction const& d) { return d; }
    };
    
      
    // GradientOf
    template<typename Id, typename Vector, typename Name, typename Direction>
    struct Diff< Id, GradientOf<Vector, Name>, Direction >
    {
      typedef GradientOf<Vector, Name> Term;
      typedef CValue<0> type;
      typedef CValue<0> dir_type;
      
      static type eval(Term const& t) { return CValue<0>(); }    
      static dir_type eval(Term const& t, Direction const& d) { return CValue<0>(); }
    };
    
    template<typename Id, typename Vector, typename Direction>
    struct Diff< Id, GradientOf<Vector, Id>, Direction >
    {
      typedef GradientOf<Vector, Id> Term;
      typedef CValue<0> type;
      typedef GradientOf<Vector, typename Direction::id> dir_type;
      
      static type eval(Term const& t) { return CValue<0>(); }    
      static dir_type eval(Term const& t, Direction const& d) { return gradientOf<typename Direction::id>(d.vecDV); }
    };
    
      
    // DerivativeOf
    template<typename Id, int I, typename Vector, typename Name, typename Direction>
    struct Diff< Id, DerivativeOf<I, Vector, Name>, Direction >
    {
      typedef DerivativeOf<I, Vector, Name> Term;
      typedef CValue<0> type;
      typedef CValue<0> dir_type;
      
      static type eval(Term const& t) { return CValue<0>(); }
      static dir_type eval(Term const& t, Direction const& d) { return CValue<0>(); }
    };
    
    template<typename Id, int I, typename Vector, typename Direction>
    struct Diff< Id, DerivativeOf<I, Vector, Id>, Direction >
    {
      typedef DerivativeOf<I, Vector, Id>                      original_type;
      typedef CValue<0>                                        type;
      typedef DerivativeOf<I, Vector, typename Direction::id>  dir_type;
      
      static type eval(original_type const& t) { return CValue<0>(); }
      static dir_type eval(original_type const& t, Direction const& d) { return derivativeOf<typename Direction::id, I>(d.vecDV); }
    };
    
  } // end namespace expressions

    
  // _____________________________________________________________________________  
  // multiplication expressions
  namespace expressions 
  {
      
    template<typename Id, typename Term1, typename Term2, typename Direction>
    class Diff< Id, Mult<Term1, Term2>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term1>::type >::type D1;
      typedef typename Simplify< typename Diff<Id, Term2>::type >::type D2;
      typedef typename Simplify< typename Diff<Id, Term1, Direction>::dir_type >::type D1_;
      typedef typename Simplify< typename Diff<Id, Term2, Direction>::dir_type >::type D2_;
      
    public:
      typedef Mult<Term1, Term2> original_type;
      typedef typename Simplify< Add< Mult< D1, Term2 >, 
				      Mult< Term1, D2 > > >::type type;
      typedef typename Simplify< Add< Mult< D1_, Term2 >, 
				      Mult< Term1, D2_> > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(simplify(diff<Id>(t.term1)) * t.term2 + t.term1 * simplify(diff<Id>(t.term2)));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term1,d)) * t.term2 + t.term1 * simplify(diff<Id>(t.term2,d)));
      }
    };
    
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, MultInverse<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef MultInverse<Term> original_type;
      typedef typename Simplify< Mult< Negative<D>, 
				      MultInverse< Pow<2, Term> > > >::type type;
      typedef typename Simplify< Mult< Negative<D_>, 
				      MultInverse< Pow<2, Term> > > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify((-simplify(diff<Id>(t.term))) / pow<2>(t.term));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify((-simplify(diff<Id>(t.term,d))) / pow<2>(t.term));
      }
    };
    
  } // end namespace expressions


  // _____________________________________________________________________________ 
  // add expressions
  namespace expressions 
  {
      
    template<typename Id, typename Term1, typename Term2, typename Direction>
    class Diff< Id, Add<Term1, Term2>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term1>::type >::type D1;
      typedef typename Simplify< typename Diff<Id, Term2>::type >::type D2;
      typedef typename Simplify< typename Diff<Id, Term1, Direction>::dir_type >::type D1_;
      typedef typename Simplify< typename Diff<Id, Term2, Direction>::dir_type >::type D2_;
      
    public:
      typedef Add<Term1, Term2> original_type;
      typedef typename Simplify< Add< D1, D2 > >::type type;
      typedef typename Simplify< Add< D1_, D2_ > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(simplify(diff<Id>(t.term1)) + simplify(diff<Id>(t.term2)));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term1,d)) + simplify(diff<Id>(t.term2,d)));
      }
    };
    
    template<typename Id, typename Term, typename Direction>
    struct Diff< Id, Negative<Term>, Direction >
    {
      typedef Negative<Term> original_type;
      typedef typename Simplify< Negative< typename Diff<Id, Term>::type > >::type type;
      typedef typename Simplify< Negative< typename Diff<Id, Term, Direction>::dir_type > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(-diff<Id>(t.term));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(-diff<Id>(t.term,d));
      }
    };
    
  } // end namespace expressions


  // _____________________________________________________________________________
  // absolute value and signum
  namespace expressions 
  {
    template<typename Term, typename Id>
    struct Diff< Id, Abs<Term> > {}; // not yet implemented
    
  } // end namespace expressions
    
    
  // _____________________________________________________________________________
  // powers and roots
  namespace expressions 
  {    
    template<typename Id, int I, typename Term, typename Direction>
    class Diff< Id, Pow<I, Term>,Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Pow<I, Term> original_type;
      typedef typename Simplify< Mult< Mult<D,  CValue<I> >, Pow<I-1, Term> > >::type type;
      typedef typename Simplify< Mult< Mult<D_,  CValue<I> >, Pow<I-1, Term> > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) * CValue<I>() * pow<I-1>(t.term));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) * CValue<I>() * pow<I-1>(t.term));
      }
    };
    
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Pow<1, Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Pow<1, Term> original_type;
      typedef D type;
      typedef D_ dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(diff<Id>(t.term));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(diff<Id>(t.term, d));
      }
    };
    
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Sqrt<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Sqrt<Term> original_type;
      typedef typename Simplify< Mult< D, MultInverse< Mult< CValue<2>, Sqrt<Term> > > > >::type type;
      typedef typename Simplify< Mult< D_, MultInverse< Mult< CValue<2>, Sqrt<Term> > > > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) / (CValue<2>() * sqrt(t.term)));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) / (CValue<2>() * sqrt(t.term)));
      }
    };
    
  } // end namespace expressions


  // _____________________________________________________________________________
  // exponential function and logarithms
  namespace expressions 
  {    
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Exp<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Exp<Term> original_type;
      typedef typename Simplify< Mult< D, Exp<Term> > >::type type;
      typedef typename Simplify< Mult< D_, Exp<Term> > >::type dir_type;
				      
      static type eval(Exp<Term> const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) * exp(t.term));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) * exp(t.term));
      }
    };
      
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Log<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Log<Term> original_type;
      typedef typename Simplify< Mult< D, MultInverse<Term> > >::type type;
      typedef typename Simplify< Mult< D_, MultInverse<Term> > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) / t.term);
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) / t.term);
      }
    };
    
  } // end namespace expressions
      
  // _____________________________________________________________________________     
  // cosine, sine and tangence
  namespace expressions 
  {    
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Cos<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Cos<Term> original_type;
      typedef typename Simplify< Negative< Mult< D, Sin<Term> > > >::type type;
      typedef typename Simplify< Negative< Mult< D_, Sin<Term> > > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(-(simplify(diff<Id>(t.term)) * sin(t.term)));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(-(simplify(diff<Id>(t.term, d)) * sin(t.term)));
      }
    };
      
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Sin<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Sin<Term> original_type;
      typedef typename Simplify< Mult< D, Cos<Term> > >::type type;
      typedef typename Simplify< Mult< D_, Cos<Term> > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) * cos(t.term));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) * cos(t.term));
      }
    };
      
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Tan<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Tan<Term> original_type;
      typedef typename Simplify< Mult< D, Add< CValue<1>, Negative< Pow<2, Tan<Term> > > > > >::type type;
      typedef typename Simplify< Mult< D_, Add< CValue<1>, Negative< Pow<2, Tan<Term> > > > > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) * (CValue<1>() - pow<2>(tan(t.term))));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) * (CValue<1>() - pow<2>(tan(t.term))));
      }
    };
    
  } // end namespace expressions

  // _____________________________________________________________________________    
  // arcus-(cosine, sine, tangence)
  namespace expressions 
  {    
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Acos<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Acos<Term> original_type;
      typedef typename Simplify< Negative< Mult< D, MultInverse< Add< CValue<1>, Negative< Pow<2, Term> > > > > > >::type type;
      typedef typename Simplify< Negative< Mult< D_, MultInverse< Add< CValue<1>, Negative< Pow<2, Term> > > > > > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(-(simplify(diff<Id>(t.term)) / (CValue<1>() - pow<2>(t.term))));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(-(simplify(diff<Id>(t.term, d)) / (CValue<1>() - pow<2>(t.term))));
      }
    };
      
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Asin<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Asin<Term> original_type;
      typedef typename Simplify< Mult< D, MultInverse< Sqrt< Add< CValue<1>, Negative< Pow<2, Term> > > > > > >::type type;
      typedef typename Simplify< Mult< D_, MultInverse< Sqrt< Add< CValue<1>, Negative< Pow<2, Term> > > > > > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) / sqrt(CValue<1>() - pow<2>(t.term)));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) / sqrt(CValue<1>() - pow<2>(t.term)));
      }
    };
      
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Atan<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Atan<Term> original_type;
      typedef typename Simplify< Mult< D, MultInverse< Add< CValue<1>, Pow<2, Term> > > > >::type type;
      typedef typename Simplify< Mult< D_, MultInverse< Add< CValue<1>, Pow<2, Term> > > > >::type dir_type;
				      
      static type eval(original_type const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) / (CValue<1>() + pow<2>(t.term)));
      }  
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) / (CValue<1>() + pow<2>(t.term)));
      }
    };
      
    template<typename Id, typename Term1, typename Term2>
    struct Diff< Id, Atan2<Term1, Term2> > {};
    
  } // end namespace expressions

  // _____________________________________________________________________________    
  // (cose, sine ,tangence)-hyperbolicus
  namespace expressions 
  {    
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Cosh<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Cosh<Term> original_type;
      typedef typename Simplify< Mult< D, Sinh<Term> > >::type type;
      typedef typename Simplify< Mult< D_, Sinh<Term> > >::type dir_type;
				      
      static type eval(Cosh<Term> const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) * sinh(t.term));
      }
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) * sinh(t.term));
      }
    };
      
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Sinh<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Sinh<Term> original_type;
      typedef typename Simplify< Mult< D, Cosh<Term> > >::type type;
      typedef typename Simplify< Mult< D_, Cosh<Term> > >::type dir_type;
				      
      static type eval(Sinh<Term> const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) * cosh(t.term));
      }
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) * cosh(t.term));
      }
    };
      
    template<typename Id, typename Term, typename Direction>
    class Diff< Id, Tanh<Term>, Direction >
    {
      typedef typename Simplify< typename Diff<Id, Term>::type >::type D;
      typedef typename Simplify< typename Diff<Id, Term, Direction>::dir_type >::type D_;
      
    public:
      typedef Tanh<Term> original_type;
      typedef typename Simplify< Mult< D, Add< CValue<1>, Negative< Pow<2, Tanh<Term> > > > > >::type type;
      typedef typename Simplify< Mult< D_, Add< CValue<1>, Negative< Pow<2, Tanh<Term> > > > > >::type dir_type;
				      
      static type eval(Tanh<Term> const& t)
      {
	return simplify(simplify(diff<Id>(t.term)) * (CValue<1>() - pow<2>(tanh(t.term))));
      }
      
      static dir_type eval(original_type const& t, Direction const& d)
      {
	return simplify(simplify(diff<Id>(t.term, d)) * (CValue<1>() - pow<2>(tanh(t.term))));
      }
    };
    
  } // end namespace expressions


  #if 0
  // _____________________________________________________________________________
  // arcus-(cosine, sine, tangence)-hyperbolicus
  namespace expressions {
    
    template<typename Id, typename Term>
    struct Diff< Id, Acosh<Term> > {}; // not yet implemented
      
    template<typename Id, typename Term>
    struct Diff< Id, Asinh<Term> > {}; // not yet implemented
      
    template<typename Id, typename Term>
    struct Diff< Id, Atanh<Term> > {}; // not yet implemented
    
  } // end namespace expressions
    
    
  // _____________________________________________________________________________
  // maximum and mininmum
  namespace expressions {
    
    template<typename Id, typename Term1, typename Term2>
    struct Diff< Id, Max<Term1, Term2> > {};  // not yet implemented
    
    template<typename Id, typename Term1, typename Term2>
    struct Diff< Id, Min<Term1, Term2> > {};  // not yet implemented
      
  } // end namespace expressions
  #endif   

  // =============================================================================
  // higher order derivatives

  // forward declaration
  namespace expressions 
  {
    template<int N, class Id, class Term> struct DiffN;
    
  } // end namespace expressions

  template<int N, typename Id, typename Term>
  inline typename expressions::DiffN<N,Id,Term>::type diff(const Term& t) { return expressions::DiffN<N,Id,Term>::eval(t); }

  namespace expressions 
  {    
    template<int N, class Id, class Term>
    struct DiffN
    {
      typedef typename DiffN<N-1, Id, typename Diff<Id, Term>::type>::type type;
			      
      static type eval(Term const& t)
      {
	return diff<N-1, Id>(diff<Id>(t));
      }
    };
    
    template<class Id, class Term>
    struct DiffN<1, Id, Term>
    {
      typedef typename Diff<Id, Term>::type type;
			      
      static type eval(Term const& t)
      {
	return diff<Id>(t);
      }
    };
    
    template<class Id, class Term>
    struct DiffN<0, Id, Term>
    {
      typedef Term type;
			      
      static type eval(Term const& t)
      {
	return t;
      }
    };
      
  } // end namespace expressions


} // end namespace AMDiS

#endif // AMDIS_DIFF_EXPRESSION_HPP
