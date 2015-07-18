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



/** \file cmath_expr.hpp */

#ifndef AMDIS_CMATH_EXPRESSION_HPP
#define AMDIS_CMATH_EXPRESSION_HPP

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"

#include "operations/functors.hpp"

namespace AMDiS 
{
  namespace expressions 
  {
    /// Expression that represents an absolute value |v|
    template<typename Term>
    struct Abs : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      
      Abs(const Term& term_) : super(term_) {}
      
      int getDegree() const
      {
	return super::term.getDegree();
      }

      inline value_type operator()(const int& iq) const { return std::abs(super::term(iq)); }
      
      std::string str() const { return std::string("abs(") + super::term.str() + ")"; }
    };
    
    
    /// Expressions that represents the sign of a value == sign(v)
    template<typename Term>
    struct Signum : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      
      Signum(const Term& term_) : super(term_) {}
      
      int getDegree() const
      {
	return 2*super::term.getDegree();
      }

      inline value_type operator()(const int& iq) const { return (super::term(iq) > 0.0 ? 1.0 : -1.0); }
      
      std::string str() const { return std::string("sign(") + super::term.str() + ")"; }
    };
    
    
    /// Expressions for ceiling a value
    template<typename Term>
    struct Ceil : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      
      Ceil(const Term& term_) : super(term_) {}
      
      int getDegree() const { return super::term.getDegree(); }

      inline value_type operator()(const int& iq) const { return std::ceil(super::term(iq)); }
      
      std::string str() const { return std::string("ceil(") + super::term.str() + ")"; }
    };
    
    
    /// Expressions to round to a lower integer
    template<typename Term>
    struct Floor : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      
      Floor(const Term& term_) : super(term_) {}
      
      int getDegree() const { return super::term.getDegree(); }

      inline value_type operator()(const int& iq) const { return std::floor(super::term(iq)); }
      
      std::string str() const { return std::string("floor(") + super::term.str() + ")"; }
    };
    
    
    // ___________________________________________________________________________
    
    
    /// Expressions that represents the Ith power of E
    template<int I, typename Term>
    struct Pow : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename AMDiS::detail::Pow<I,typename Term::value_type>::result_type value_type;
      
      Pow(const Term& term_) : super(term_) {}
      
      int getDegree() const
      {
	return I * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return functors::pow<I,typename Term::value_type>::eval(super::term(iq)); }
      
      std::string str() const { return std::string("pow<") + boost::lexical_cast<std::string>(I) + ">(" + super::term.str() + ")"; }
    };

    
    /// Expression that represents that square root of E
    template<typename Term>
    struct Sqrt : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      
      Sqrt(const Term& term_) : super(term_) {}
      
      int getDegree() const
      {
	return 2 * (super::term.getDegree()); // stimmt nicht ganz
      }

      inline value_type operator()(const int& iq) const { return std::sqrt(super::term(iq)); }
      
      std::string str() const { return std::string("sqrt(") + super::term.str() + ")"; }
    };

    
    // ___________________________________________________________________________  
    
    /// Expressions for the exponential function
    template<typename Term>
    struct Exp : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Exp(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::exp(super::term(iq)); }
      
      std::string str() const { return std::string("exp(") + super::term.str() + ")"; }
    };

    /// Expression for the logarithm
    template<typename Term>
    struct Log : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Log(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::log(super::term(iq)); }
      
      std::string str() const { return std::string("log(") + super::term.str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for Cosine
    template<typename Term>
    struct Cos : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Cos(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::cos(super::term(iq)); }
      
      std::string str() const { return std::string("cos(") + super::term.str() + ")"; }
    };

    /// Expressions for Sine
    template<typename Term>
    struct Sin : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Sin(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::sin(super::term(iq)); }
      
      std::string str() const { return std::string("sin(") + super::term.str() + ")"; }
    };

    /// Expressions for Tangence
    template<typename Term>
    struct Tan : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Tan(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::tan(super::term(iq)); }
      
      std::string str() const { return std::string("tan(") + super::term.str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for arcus cosine
    template<typename Term>
    struct Acos : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Acos(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::acos(super::term(iq)); }
      
      std::string str() const { return std::string("acos(") + super::term.str() + ")"; }
    };

    /// Expressions for arcus sine
    template<typename Term>
    struct Asin : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Asin(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::asin(super::term(iq)); }
      
      std::string str() const { return std::string("asin(") + super::term.str() + ")"; }
    };

    /// Expressions for arcus tangence
    template<typename Term>
    struct Atan : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Atan(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::atan(super::term(iq)); }
      
      std::string str() const { return std::string("atan(") + super::term.str() + ")"; }
    };

    /// Expressions for arcus tangence2, i.e. atan(x/y)
    template<typename Term1, typename Term2>
    struct Atan2 : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> super;
      typedef typename Term1::value_type value_type;
      int degree;
      
      Atan2(const Term1& term1_, const Term2& term2_, int degree_ = 1) : super(term1_, term2_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term1.getDegree() + super::term2.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::atan2(super::term1(iq), super::term2(iq)); }
      
      std::string str() const { return std::string("atan2(") + super::term.str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for cosine hyperbolicus
    template<typename Term>
    struct Cosh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Cosh(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::cosh(super::term(iq)); }
      
      std::string str() const { return std::string("cosh(") + super::term.str() + ")"; }
    };

    /// Expressions for sine hyperbolicus
    template<typename Term>
    struct Sinh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Sinh(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::sinh(super::term(iq)); }
      
      std::string str() const { return std::string("sinh(") + super::term.str() + ")"; }
    };

    /// Expressions for tangence hyperbolicus
    template<typename Term>
    struct Tanh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Tanh(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return std::tanh(super::term(iq)); }
      
      std::string str() const { return std::string("tanh(") + super::term.str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for arcus cosine hyperbolicus
    template<typename Term>
    struct Acosh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Acosh(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { 
	value_type tmp = super::term(iq);
	return std::log(tmp + std::sqrt(sqr(tmp) - 1.0)); 
      }
      
      std::string str() const { return std::string("acosh(") + super::term.str() + ")"; }
    };

    /// Expressions for arcus sine hyperbolicus
    template<typename Term>
    struct Asinh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Asinh(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { 
	value_type tmp = super::term(iq);
	return std::log(tmp + std::sqrt(sqr(tmp) + 1.0)); 
      }
      
      std::string str() const { return std::string("asinh(") + super::term.str() + ")"; }
    };

    /// Expressions for arcus tangence hyperbolicus
    template<typename Term>
    struct Atanh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Atanh(const Term& term_, int degree_ = 1) : super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { 
	value_type tmp = super::term(iq);
	return 0.5 * std::log((1.0 + tmp) / (1.0 - tmp)); 
      }
      
      std::string str() const { return std::string("atanh(") + super::term.str() + ")"; }
    };
    
    
    // ___________________________________________________________________________
    
    
    /// Expressions for the maximum of two expressions max(E1, E2)
    template<typename Term1, typename Term2>
    struct Max : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> super;
      typedef typename Term1::value_type value_type;
      
      Max(const Term1& term1_, const Term2& term2_)
	: super(term1_, term2_) {}
      
      int getDegree() const
      {
	return std::max(super::term1.getDegree(), super::term2.getDegree());
      }

      inline value_type operator()(const int& iq) const
      {
	return std::max(super::term1(iq), super::term2(iq)); 
      }
      
      std::string str() const { return std::string("max(") + super::term1.str() + ", " + super::term2.str() + ")"; }
    };
    
    
    /// Expressions for the minimum of two expressions min(E1, E2)
    template<typename Term1, typename Term2>
    struct Min : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> super;
      typedef typename Term1::value_type value_type;
      
      Min(const Term1& term1_, const Term2& term2_)
	: super(term1_, term2_) {}
      
      int getDegree() const
      {
	return std::max(super::term1.getDegree(), super::term2.getDegree());
      }

      inline value_type operator()(const int& iq) const 
      {
	return std::min(super::term1(iq), super::term2(iq)); 
      }
      
      std::string str() const { return std::string("min(") + super::term1.str() + ", " + super::term2.str() + ")"; }
    };
    
  } // end namespace expressions


  namespace result_of
  {
    template<typename Term1, typename Term2>
    struct Min : boost::enable_if
      < 
	typename traits::is_valid_arg2<Term1, Term2>::type,
	expressions::Min
	<
	  typename traits::to_expr<Term1>::type, 
	  typename traits::to_expr<Term2>::type
	>
      > {};
      
      
    template<typename Term1, typename Term2>
    struct Max : boost::enable_if
      < 
	typename traits::is_valid_arg2<Term1, Term2>::type,
	expressions::Max
	<
	  typename traits::to_expr<Term1>::type, 
	  typename traits::to_expr<Term2>::type
	>
      > {};
      
  } // end namespace result_of


  // maximum of two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Max<Term1, Term2>::type
  max(const Term1& t1, const Term2& t2) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    return expressions::Max< typename Expr1::type, typename Expr2::type >
	    (Expr1::get(t1), Expr2::get(t2));
  }

  // minimum of two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Min<Term1, Term2>::type
  min(const Term1& t1, const Term2& t2) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    return expressions::Min< typename Expr1::type, typename Expr2::type >
	    (Expr1::get(t1), Expr2::get(t2));
  }


  //______________________________________________________________________________

  // absolute value of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Abs<Term> >::type
  abs_(const Term& t) { return expressions::Abs<Term>(t); } // TODO: Funktionsnamen ohne Unterstrich

  // signum of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Signum<Term> >::type
  signum(const Term& t) { return expressions::Signum<Term>(t); }

  // 
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Ceil<Term> >::type
  ceil(const Term& t) { return expressions::Ceil<Term>(t); }

  // 
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Floor<Term> >::type
  floor(const Term& t) { return expressions::Floor<Term>(t); }

  //______________________________________________________________________________

  // I'th power of a term
  template<int I, typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Pow<I, Term> >::type
  pow(const Term& t) { return expressions::Pow<I, Term>(t); }

  // square root of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Sqrt<Term> >::type
  sqrt(const Term& t) { return expressions::Sqrt<Term>(t); }

  //______________________________________________________________________________

  // exponential function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Exp<Term> >::type
  exp(const Term& t) { return expressions::Exp<Term>(t); }

  // natural logarithm of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Log<Term> >::type
  log(const Term& t) { return expressions::Log<Term>(t); }

  //______________________________________________________________________________

  // cosine function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Cos<Term> >::type
  cos(const Term& t) { return expressions::Cos<Term>(t); }

  // sine function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Sin<Term> >::type
  sin(const Term& t) { return expressions::Sin<Term>(t); }

  // tangens function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Tan<Term> >::type
  tan(const Term& t) { return expressions::Tan<Term>(t); }

  //______________________________________________________________________________

  // arkuscosine function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Acos<Term> >::type
  acos(const Term& t) { return expressions::Acos<Term>(t); }

  // arkussine function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Asin<Term> >::type
  asin(const Term& t) { return expressions::Asin<Term>(t); }

  // arkustangens function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Atan<Term> >::type
  atan(const Term& t) { return expressions::Atan<Term>(t); }

  //______________________________________________________________________________

  // cosine-hyperbolicus function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Cosh<Term> >::type
  cosh(const Term& t) { return expressions::Cosh<Term>(t); }

  // sine-hyperbolicus function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Sinh<Term> >::type
  sinh(const Term& t) { return expressions::Sinh<Term>(t); }

  // tangens-hyperbolicus function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Tanh<Term> >::type
  tanh(const Term& t) { return expressions::Tanh<Term>(t); }

  //______________________________________________________________________________

  // arkuscosine-hypüerbolicus function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Acosh<Term> >::type
  acosh(const Term& t) { return expressions::Acosh<Term>(t); }

  // arkussine-hyperbolicus function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Asinh<Term> >::type
  asinh(const Term& t) { return expressions::Asinh<Term>(t); }

  // arkustangens-hyperbolicus function of a term
  template<typename Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Atanh<Term> >::type
  atanh(const Term& t) { return expressions::Atanh<Term>(t); }

} // end namespace AMDiS

#endif // AMDIS_CMATH_EXPRESSION_HPP
