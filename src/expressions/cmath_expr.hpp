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
    struct Abs : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      template <class Term_>
      Abs(Term_&& term_) 
        : Super(std::forward<Term_>(term_)) {}
      
      int getDegree() const
      {
        return Super::getDegree(_0);
      }

      value_type operator()(const int& iq) const { return std::abs(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("abs(") + Super::getTerm(_0).str() + ")"; }
    };
    
    
    /// Expressions that represents the sign of a value == sign(v)
    template <class Term>
    struct Signum : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      template <class Term_>
      Signum(Term_&& term_) 
        : Super(std::forward<Term_>(term_)) {}
      
      int getDegree() const
      {
        return 2*Super::getDegree(_0);
      }

      value_type operator()(const int& iq) const { return (Super::getTerm(_0)(iq) > 0.0 ? 1.0 : -1.0); }
      
      std::string str() const { return std::string("sign(") + Super::getTerm(_0).str() + ")"; }
    };
    
    
    /// Expressions for ceiling a value
    template <class Term>
    struct Ceil : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      template <class Term_>
      Ceil(Term_&& term_) 
        : Super(std::forward<Term_>(term_)) {}
      
      int getDegree() const { return Super::getDegree(_0); }

      value_type operator()(const int& iq) const { return std::ceil(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("ceil(") + Super::getTerm(_0).str() + ")"; }
    };
    
    
    /// Expressions to round to a lower integer
    template <class Term>
    struct Floor : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      template <class Term_>
      Floor(Term_&& term_) 
        : Super(std::forward<Term_>(term_)) {}
      
      int getDegree() const { return Super::getDegree(_0); }

      value_type operator()(const int& iq) const { return std::floor(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("floor(") + Super::getTerm(_0).str() + ")"; }
    };
    
    
    // ___________________________________________________________________________
    
    
    /// Expressions that represents the Ith power of E
    template <int I, class Term>
    struct Pow : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
//      typedef typename AMDiS::detail::Pow<I,typename Term::value_type>::result_type value_type;
      using value_type = Value_t<Term>;
      
      template <class Term_>
      Pow(Term_&& term_) 
        : Super(std::forward<Term_>(term_)) {}
      
      int getDegree() const
      {
        return I * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return functors::pow<I,Value_t<Term> >::eval(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("pow<") + std::to_string(I) + ">(" + Super::getTerm(_0).str() + ")"; }
    };

    
    /// Expression that represents that square root of E
    template <class Term>
    struct Sqrt : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      template <class Term_>
      Sqrt(Term_&& term_) 
        : Super(std::forward<Term_>(term_)) {}
      
      int getDegree() const
      {
        return 2 * (Super::getDegree(_0)); // stimmt nicht ganz
      }

      value_type operator()(const int& iq) const { return std::sqrt(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("sqrt(") + Super::getTerm(_0).str() + ")"; }
    };

    
    // ___________________________________________________________________________  
    
    /// Expressions for the exponential function
    template <class Term>
    struct Exp : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Exp(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::exp(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("exp(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expression for the logarithm
    template <class Term>
    struct Log : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Log(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::log(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("log(") + Super::getTerm(_0).str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for Cosine
    template <class Term>
    struct Cos : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Cos(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::cos(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("cos(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expressions for Sine
    template <class Term>
    struct Sin : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Sin(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::sin(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("sin(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expressions for Tangence
    template <class Term>
    struct Tan : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Tan(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::tan(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("tan(") + Super::getTerm(_0).str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for arcus cosine
    template <class Term>
    struct Acos : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Acos(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::acos(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("acos(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expressions for arcus sine
    template <class Term>
    struct Asin : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Asin(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::asin(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("asin(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expressions for arcus tangence
    template <class Term>
    struct Atan : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Atan(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::atan(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("atan(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expressions for arcus tangence2, i.e. atan(x/y)
    template <class Term1, class Term2>
    struct Atan2 : public LazyOperatorTerms<Term1, Term2>
    {
      using Super = LazyOperatorTerms<Term1, Term2>;
      using value_type = Value_t<Term1>;
      int deg;
      
      template <class Term0_, class Term1_>
      Atan2(Term0_&& term0_, Term1_&& term1_, int degree_ = 1) 
        : Super(std::forward<Term0_>(term0_), std::forward<Term1_>(term1_)), deg(degree_) 
      { }
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0) + Super::getDegree(_1));
      }

      value_type operator()(const int& iq) const { return std::atan2(Super::getTerm(_0)(iq), Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("atan2(") + Super::getTerm(_0).str() + ", " + Super::getTerm(_1) + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for cosine hyperbolicus
    template <class Term>
    struct Cosh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Cosh(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::cosh(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("cosh(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expressions for sine hyperbolicus
    template <class Term>
    struct Sinh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Sinh(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::sinh(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("sinh(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expressions for tangence hyperbolicus
    template <class Term>
    struct Tanh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Tanh(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (Super::getDegree(_0));
      }

      value_type operator()(const int& iq) const { return std::tanh(Super::getTerm(_0)(iq)); }
      
      std::string str() const { return std::string("tanh(") + Super::getTerm(_0).str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for arcus cosine hyperbolicus
    template <class Term>
    struct Acosh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Acosh(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * Super::getDegree(_0);
      }

      value_type operator()(const int& iq) const 
      { 
        value_type tmp = Super::getTerm(_0)(iq);
      	return std::log(tmp + std::sqrt(sqr(tmp) - 1.0)); 
      }
      
      std::string str() const { return std::string("acosh(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expressions for arcus sine hyperbolicus
    template <class Term>
    struct Asinh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Asinh(Term_&& term_, int degree_ = 1) 
        : Super(std::forward<Term_>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * Super::getDegree(_0);
      }

      value_type operator()(const int& iq) const 
      { 
        value_type tmp = Super::getTerm(_0)(iq);
        return std::log(tmp + std::sqrt(sqr(tmp) + 1.0)); 
      }
      
      std::string str() const { return std::string("asinh(") + Super::getTerm(_0).str() + ")"; }
    };

    /// Expressions for arcus tangence hyperbolicus
    template <class Term>
    struct Atanh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      template <class Term_>
      Atanh(Term_&& term_, int deg_ = 1) : Super(std::forward<Term_>(term_)), deg(deg_) {}
      
      int getDegree() const
      {
        return deg * Super::getDegree(_0);
      }

      value_type operator()(const int& iq) const 
      { 
      	value_type tmp = Super::getTerm(_0)(iq);
      	return 0.5 * std::log((1.0 + tmp) / (1.0 - tmp)); 
      }
      
      std::string str() const { return std::string("atanh(") + Super::getTerm(_0).str() + ")"; }
    };
    
    
    // ___________________________________________________________________________
    
    
    /// Expressions for the maximum of two expressions max(E1, E2)
    template <class Term1, class Term2>
    struct Max : public LazyOperatorTerms<Term1, Term2>
    {
      using Super = LazyOperatorTerms<Term1, Term2>;
      using value_type = typename std::common_type<Value_t<Term1>, Value_t<Term2> >::type;
      
      template <class Term0_, class Term1_>
      Max(Term0_&& term0_, Term1_&& term1_)
        : Super(std::forward<Term0_>(term0_), std::forward<Term1_>(term1_)) 
      { }
      
      int getDegree() const
      {
        return std::max(Super::getDegree(_0), Super::getDegree(_1));
      }

      value_type operator()(const int& iq) const
      {
        return std::max(Super::getTerm(_0)(iq), Super::getTerm(_1)(iq)); 
      }
      
      std::string str() const { return std::string("max(") + Super::getTerm(_0).str() + ", " + Super::getTerm(_1).str() + ")"; }
    };
    
    
    /// Expressions for the minimum of two expressions min(E1, E2)
    template <class Term1, class Term2>
    struct Min : public LazyOperatorTerms<Term1, Term2>
    {
      using Super = LazyOperatorTerms<Term1, Term2>;
      using value_type = typename std::common_type<Value_t<Term1>, Value_t<Term2> >::type;
      
      template <class Term0_, class Term1_>
      Min(Term0_&& term0_, Term1_&& term1_)
        : Super(std::forward<Term0_>(term0_), std::forward<Term1_>(term1_)) 
      { }
      
      int getDegree() const
      {
        return std::max(Super::getDegree(_0), Super::getDegree(_1));
      }

      value_type operator()(const int& iq) const 
      {
      	return std::min(Super::getTerm(_0)(iq), Super::getTerm(_1)(iq)); 
      }
      
      std::string str() const { return std::string("min(") + Super::getTerm(_0).str() + ", " + Super::getTerm(_1).str() + ")"; }
    };
    
  } // end namespace expressions


  namespace result_of
  {
    template <class Term1, class Term2>
    using Min = enable_if
      < 
      	traits::is_valid_arg2<Term1, Term2>,
      	expressions::Min
      	<
      	  typename traits::to_expr<Term1>::type, 
      	  typename traits::to_expr<Term2>::type
      	>
      >;
      
      
    template <class Term1, class Term2>
    using Max = enable_if
      < 
      	traits::is_valid_arg2<Term1, Term2>,
      	expressions::Max
      	<
      	  typename traits::to_expr<Term1>::type, 
      	  typename traits::to_expr<Term2>::type
      	>
      >;
      
  } // end namespace result_of


  // maximum of two terms
  // _____________________________________________________________________________
  template <class Term0, class Term1>
  inline typename result_of::Max<Term0, Term1>::type
  max(Term0&& t0, Term1&& t1) 
  { 
    using Expr0 = traits::to_expr<Term0>;
    using Expr1 = traits::to_expr<Term1>;
    return {Expr0::get(t0), Expr1::get(t1)};
  }

  // minimum of two terms
  // _____________________________________________________________________________
  template <class Term0, class Term1>
  inline typename result_of::Min<Term0, Term1>::type
  min(Term0&& t0, Term1&& t1) 
  { 
    using Expr0 = traits::to_expr<Term0>;
    using Expr1 = traits::to_expr<Term1>;
    return {Expr0::get(t0), Expr1::get(t1)};
  }


  //______________________________________________________________________________

  // absolute value of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Abs<Term> >::type
  abs_(Term&& t) { return {std::forward<Term>(t)}; } // TODO: Funktionsnamen ohne Unterstrich

  // signum of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Signum<Term> >::type
  signum(Term&& t) { return {std::forward<Term>(t)}; }

  // 
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Ceil<Term> >::type
  ceil(Term&& t) { return {std::forward<Term>(t)}; }

  // 
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Floor<Term> >::type
  floor(Term&& t) { return {std::forward<Term>(t)}; }

  //______________________________________________________________________________

  // I'th power of a term
  template <int I, class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Pow<I, Term> >::type
  pow(Term&& t) { return {std::forward<Term>(t)}; }

  // square root of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Sqrt<Term> >::type
  sqrt(Term&& t) { return {std::forward<Term>(t)}; }

  //______________________________________________________________________________

  // exponential function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Exp<Term> >::type
  exp(Term&& t) { return {std::forward<Term>(t)}; }

  // natural logarithm of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Log<Term> >::type
  log(Term&& t) { return {std::forward<Term>(t)}; }

  //______________________________________________________________________________

  // cosine function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Cos<Term> >::type
  cos(Term&& t) { return {std::forward<Term>(t)}; }

  // sine function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Sin<Term> >::type
  sin(Term&& t) { return {std::forward<Term>(t)}; }

  // tangens function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Tan<Term> >::type
  tan(Term&& t) { return {std::forward<Term>(t)}; }

  //______________________________________________________________________________

  // arkuscosine function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Acos<Term> >::type
  acos(Term&& t) { return {std::forward<Term>(t)}; }

  // arkussine function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Asin<Term> >::type
  asin(Term&& t) { return {std::forward<Term>(t)}; }

  // arkustangens function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Atan<Term> >::type
  atan(Term&& t) { return {std::forward<Term>(t)}; }

  //______________________________________________________________________________

  // cosine-hyperbolicus function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Cosh<Term> >::type
  cosh(Term&& t) { return {std::forward<Term>(t)}; }

  // sine-hyperbolicus function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Sinh<Term> >::type
  sinh(Term&& t) { return {std::forward<Term>(t)}; }

  // tangens-hyperbolicus function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Tanh<Term> >::type
  tanh(Term&& t) { return {std::forward<Term>(t)}; }

  //______________________________________________________________________________

  // arkuscosine-hypüerbolicus function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Acosh<Term> >::type
  acosh(Term&& t) { return {std::forward<Term>(t)}; }

  // arkussine-hyperbolicus function of a term
  template <class Term>
  inline typename enable_if< traits::is_expr<Term>,
    expressions::Asinh<Term> >::type
  asinh(Term&& t) { return {std::forward<Term>(t)}; }

  // arkustangens-hyperbolicus function of a term
  template<class Term>
  inline typename enable_if<
    traits::is_expr<Term>,
    expressions::Atanh<Term> >::type
  atanh(Term&& t) { return {std::forward<Term>(t)}; }

} // end namespace AMDiS
