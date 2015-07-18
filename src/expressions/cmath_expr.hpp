/** \file cmath_expr.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"

#include "operations/functors.hpp"

namespace AMDiS 
{
  namespace expressions 
  {
    /// Expression that represents an absolute value |v|
    template <class Term>
    struct Abs : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      
      Abs(const Term& term_) : Super(term_) {}
      
      int getDegree() const
      {
	return Super::term.getDegree();
      }

      value_type operator()(const int& iq) const { return std::abs(Super::term(iq)); }
      
      std::string str() const { return std::string("abs(") + Super::term.str() + ")"; }
    };
    
    
    /// Expressions that represents the sign of a value == sign(v)
    template <class Term>
    struct Signum : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      
      Signum(const Term& term_) : Super(term_) {}
      
      int getDegree() const
      {
	return 2*Super::term.getDegree();
      }

      value_type operator()(const int& iq) const { return (Super::term(iq) > 0.0 ? 1.0 : -1.0); }
      
      std::string str() const { return std::string("sign(") + Super::term.str() + ")"; }
    };
    
    
    /// Expressions for ceiling a value
    template <class Term>
    struct Ceil : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      
      Ceil(const Term& term_) : Super(term_) {}
      
      int getDegree() const { return Super::term.getDegree(); }

      value_type operator()(const int& iq) const { return std::ceil(Super::term(iq)); }
      
      std::string str() const { return std::string("ceil(") + Super::term.str() + ")"; }
    };
    
    
    /// Expressions to round to a lower integer
    template <class Term>
    struct Floor : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      
      Floor(const Term& term_) : Super(term_) {}
      
      int getDegree() const { return Super::term.getDegree(); }

      value_type operator()(const int& iq) const { return std::floor(Super::term(iq)); }
      
      std::string str() const { return std::string("floor(") + Super::term.str() + ")"; }
    };
    
    
    // ___________________________________________________________________________
    
    
    /// Expressions that represents the Ith power of E
    template <int I, class Term>
    struct Pow : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename AMDiS::detail::Pow<I,typename Term::value_type>::result_type value_type;
      
      Pow(const Term& term_) : Super(term_) {}
      
      int getDegree() const
      {
	return I * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return functors::pow<I,typename Term::value_type>::eval(Super::term(iq)); }
      
      std::string str() const { return std::string("pow<") + std::to_string(I) + ">(" + Super::term.str() + ")"; }
    };

    
    /// Expression that represents that square root of E
    template <class Term>
    struct Sqrt : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      
      Sqrt(const Term& term_) : Super(term_) {}
      
      int getDegree() const
      {
	return 2 * (Super::term.getDegree()); // stimmt nicht ganz
      }

      value_type operator()(const int& iq) const { return std::sqrt(Super::term(iq)); }
      
      std::string str() const { return std::string("sqrt(") + Super::term.str() + ")"; }
    };

    
    // ___________________________________________________________________________  
    
    /// Expressions for the exponential function
    template <class Term>
    struct Exp : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Exp(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::exp(Super::term(iq)); }
      
      std::string str() const { return std::string("exp(") + Super::term.str() + ")"; }
    };

    /// Expression for the logarithm
    template <class Term>
    struct Log : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Log(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::log(Super::term(iq)); }
      
      std::string str() const { return std::string("log(") + Super::term.str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for Cosine
    template <class Term>
    struct Cos : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Cos(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::cos(Super::term(iq)); }
      
      std::string str() const { return std::string("cos(") + Super::term.str() + ")"; }
    };

    /// Expressions for Sine
    template <class Term>
    struct Sin : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Sin(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::sin(Super::term(iq)); }
      
      std::string str() const { return std::string("sin(") + Super::term.str() + ")"; }
    };

    /// Expressions for Tangence
    template <class Term>
    struct Tan : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Tan(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::tan(Super::term(iq)); }
      
      std::string str() const { return std::string("tan(") + Super::term.str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for arcus cosine
    template <class Term>
    struct Acos : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Acos(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::acos(Super::term(iq)); }
      
      std::string str() const { return std::string("acos(") + Super::term.str() + ")"; }
    };

    /// Expressions for arcus sine
    template <class Term>
    struct Asin : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Asin(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::asin(Super::term(iq)); }
      
      std::string str() const { return std::string("asin(") + Super::term.str() + ")"; }
    };

    /// Expressions for arcus tangence
    template <class Term>
    struct Atan : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Atan(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::atan(Super::term(iq)); }
      
      std::string str() const { return std::string("atan(") + Super::term.str() + ")"; }
    };

    /// Expressions for arcus tangence2, i.e. atan(x/y)
    template <class Term1, class Term2>
    struct Atan2 : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> Super;
      typedef typename Term1::value_type value_type;
      int degree;
      
      Atan2(const Term1& term1_, const Term2& term2_, int degree_ = 1) : Super(term1_, term2_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term1.getDegree() + Super::term2.getDegree());
      }

      value_type operator()(const int& iq) const { return std::atan2(Super::term1(iq), Super::term2(iq)); }
      
      std::string str() const { return std::string("atan2(") + Super::term.str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for cosine hyperbolicus
    template <class Term>
    struct Cosh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Cosh(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::cosh(Super::term(iq)); }
      
      std::string str() const { return std::string("cosh(") + Super::term.str() + ")"; }
    };

    /// Expressions for sine hyperbolicus
    template <class Term>
    struct Sinh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Sinh(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::sinh(Super::term(iq)); }
      
      std::string str() const { return std::string("sinh(") + Super::term.str() + ")"; }
    };

    /// Expressions for tangence hyperbolicus
    template <class Term>
    struct Tanh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Tanh(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const { return std::tanh(Super::term(iq)); }
      
      std::string str() const { return std::string("tanh(") + Super::term.str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for arcus cosine hyperbolicus
    template <class Term>
    struct Acosh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Acosh(const Term& term_, int degree_ = 1) 
	: Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const 
      { 
	value_type tmp = Super::term(iq);
	return std::log(tmp + std::sqrt(sqr(tmp) - 1.0)); 
      }
      
      std::string str() const { return std::string("acosh(") + Super::term.str() + ")"; }
    };

    /// Expressions for arcus sine hyperbolicus
    template <class Term>
    struct Asinh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Asinh(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const 
      { 
	value_type tmp = Super::term(iq);
	return std::log(tmp + std::sqrt(sqr(tmp) + 1.0)); 
      }
      
      std::string str() const { return std::string("asinh(") + Super::term.str() + ")"; }
    };

    /// Expressions for arcus tangence hyperbolicus
    template <class Term>
    struct Atanh : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> Super;
      typedef typename Term::value_type value_type;
      int degree;
      
      Atanh(const Term& term_, int degree_ = 1) : Super(term_), degree(degree_) {}
      
      int getDegree() const
      {
	return degree * (Super::term.getDegree());
      }

      value_type operator()(const int& iq) const 
      { 
	value_type tmp = Super::term(iq);
	return 0.5 * std::log((1.0 + tmp) / (1.0 - tmp)); 
      }
      
      std::string str() const { return std::string("atanh(") + Super::term.str() + ")"; }
    };
    
    
    // ___________________________________________________________________________
    
    
    /// Expressions for the maximum of two expressions max(E1, E2)
    template <class Term1, class Term2>
    struct Max : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> Super;
      typedef typename Term1::value_type value_type;
      
      Max(const Term1& term1_, const Term2& term2_)
	: Super(term1_, term2_) {}
      
      int getDegree() const
      {
	return std::max(Super::term1.getDegree(), Super::term2.getDegree());
      }

      value_type operator()(const int& iq) const
      {
	return std::max(Super::term1(iq), Super::term2(iq)); 
      }
      
      std::string str() const { return std::string("max(") + Super::term1.str() + ", " + Super::term2.str() + ")"; }
    };
    
    
    /// Expressions for the minimum of two expressions min(E1, E2)
    template <class Term1, class Term2>
    struct Min : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> Super;
      typedef typename Term1::value_type value_type;
      
      Min(const Term1& term1_, const Term2& term2_)
	: Super(term1_, term2_) {}
      
      int getDegree() const
      {
	return std::max(Super::term1.getDegree(), Super::term2.getDegree());
      }

      value_type operator()(const int& iq) const 
      {
	return std::min(Super::term1(iq), Super::term2(iq)); 
      }
      
      std::string str() const { return std::string("min(") + Super::term1.str() + ", " + Super::term2.str() + ")"; }
    };
    
  } // end namespace expressions


  namespace result_of
  {
    template <class Term1, class Term2>
    struct Min : boost::enable_if
      < 
	traits::is_valid_arg2<Term1, Term2>,
	expressions::Min
	<
	  typename traits::to_expr<Term1>::type, 
	  typename traits::to_expr<Term2>::type
	>
      > {};
      
      
    template <class Term1, class Term2>
    struct Max : boost::enable_if
      < 
	traits::is_valid_arg2<Term1, Term2>,
	expressions::Max
	<
	  typename traits::to_expr<Term1>::type, 
	  typename traits::to_expr<Term2>::type
	>
      > {};
      
  } // end namespace result_of


  // maximum of two terms
  // _____________________________________________________________________________
  template <class Term1, class Term2>
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
  template <class Term1, class Term2>
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
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Abs<Term> >::type
  abs_(const Term& t) { return expressions::Abs<Term>(t); } // TODO: Funktionsnamen ohne Unterstrich

  // signum of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Signum<Term> >::type
  signum(const Term& t) { return expressions::Signum<Term>(t); }

  // 
  template <class Term>
  inline typename boost::enable_if<
    typename traits::is_expr<Term>::type,
    expressions::Ceil<Term> >::type
  ceil(const Term& t) { return expressions::Ceil<Term>(t); }

  // 
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Floor<Term> >::type
  floor(const Term& t) { return expressions::Floor<Term>(t); }

  //______________________________________________________________________________

  // I'th power of a term
  template <int I, class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Pow<I, Term> >::type
  pow(const Term& t) { return expressions::Pow<I, Term>(t); }

  // square root of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Sqrt<Term> >::type
  sqrt(const Term& t) { return expressions::Sqrt<Term>(t); }

  //______________________________________________________________________________

  // exponential function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Exp<Term> >::type
  exp(const Term& t) { return expressions::Exp<Term>(t); }

  // natural logarithm of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Log<Term> >::type
  log(const Term& t) { return expressions::Log<Term>(t); }

  //______________________________________________________________________________

  // cosine function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Cos<Term> >::type
  cos(const Term& t) { return expressions::Cos<Term>(t); }

  // sine function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Sin<Term> >::type
  sin(const Term& t) { return expressions::Sin<Term>(t); }

  // tangens function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Tan<Term> >::type
  tan(const Term& t) { return expressions::Tan<Term>(t); }

  //______________________________________________________________________________

  // arkuscosine function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Acos<Term> >::type
  acos(const Term& t) { return expressions::Acos<Term>(t); }

  // arkussine function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Asin<Term> >::type
  asin(const Term& t) { return expressions::Asin<Term>(t); }

  // arkustangens function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Atan<Term> >::type
  atan(const Term& t) { return expressions::Atan<Term>(t); }

  //______________________________________________________________________________

  // cosine-hyperbolicus function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Cosh<Term> >::type
  cosh(const Term& t) { return expressions::Cosh<Term>(t); }

  // sine-hyperbolicus function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Sinh<Term> >::type
  sinh(const Term& t) { return expressions::Sinh<Term>(t); }

  // tangens-hyperbolicus function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Tanh<Term> >::type
  tanh(const Term& t) { return expressions::Tanh<Term>(t); }

  //______________________________________________________________________________

  // arkuscosine-hypüerbolicus function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Acosh<Term> >::type
  acosh(const Term& t) { return expressions::Acosh<Term>(t); }

  // arkussine-hyperbolicus function of a term
  template <class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Asinh<Term> >::type
  asinh(const Term& t) { return expressions::Asinh<Term>(t); }

  // arkustangens-hyperbolicus function of a term
  template<class Term>
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Atanh<Term> >::type
  atanh(const Term& t) { return expressions::Atanh<Term>(t); }

} // end namespace AMDiS
