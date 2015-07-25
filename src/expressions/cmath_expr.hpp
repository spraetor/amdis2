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
      
      Abs(Term&& term_) : Super(std::forward<Term>(term_)) {}
      
      int getDegree() const
      {
        return degree<0>(*this);
      }

      value_type operator()(const int& iq) const { return std::abs(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("abs(") + term<0>(*this).str() + ")"; }
    };
    
    
    /// Expressions that represents the sign of a value == sign(v)
    template <class Term>
    struct Signum : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      Signum(Term&& term_) : Super(std::forward<Term>(term_)) {}
      
      int getDegree() const
      {
        return 2*degree<0>(*this);
      }

      value_type operator()(const int& iq) const { return (term<0>(*this)(iq) > 0.0 ? 1.0 : -1.0); }
      
      std::string str() const { return std::string("sign(") + term<0>(*this).str() + ")"; }
    };
    
    
    /// Expressions for ceiling a value
    template <class Term>
    struct Ceil : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      Ceil(Term&& term_) : Super(std::forward<Term>(term_)) {}
      
      int getDegree() const { return degree<0>(*this); }

      value_type operator()(const int& iq) const { return std::ceil(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("ceil(") + term<0>(*this).str() + ")"; }
    };
    
    
    /// Expressions to round to a lower integer
    template <class Term>
    struct Floor : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      Floor(Term&& term_) : Super(std::forward<Term>(term_)) {}
      
      int getDegree() const { return degree<0>(*this); }

      value_type operator()(const int& iq) const { return std::floor(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("floor(") + term<0>(*this).str() + ")"; }
    };
    
    
    // ___________________________________________________________________________
    
    
    /// Expressions that represents the Ith power of E
    template <int I, class Term>
    struct Pow : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
//      typedef typename AMDiS::detail::Pow<I,typename Term::value_type>::result_type value_type;
      using value_type = Value_t<Term>;
      
      Pow(Term&& term_) : Super(std::forward<Term>(term_)) {}
      
      int getDegree() const
      {
        return I * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return functors::pow<I,Value_t<Term> >::eval(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("pow<") + std::to_string(I) + ">(" + term<0>(*this).str() + ")"; }
    };

    
    /// Expression that represents that square root of E
    template <class Term>
    struct Sqrt : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      Sqrt(Term&& term_) : Super(std::forward<Term>(term_)) {}
      
      int getDegree() const
      {
        return 2 * (degree<0>(*this)); // stimmt nicht ganz
      }

      value_type operator()(const int& iq) const { return std::sqrt(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("sqrt(") + term<0>(*this).str() + ")"; }
    };

    
    // ___________________________________________________________________________  
    
    /// Expressions for the exponential function
    template <class Term>
    struct Exp : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Exp(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::exp(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("exp(") + term<0>(*this).str() + ")"; }
    };

    /// Expression for the logarithm
    template <class Term>
    struct Log : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Log(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::log(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("log(") + term<0>(*this).str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for Cosine
    template <class Term>
    struct Cos : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Cos(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::cos(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("cos(") + term<0>(*this).str() + ")"; }
    };

    /// Expressions for Sine
    template <class Term>
    struct Sin : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Sin(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::sin(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("sin(") + term<0>(*this).str() + ")"; }
    };

    /// Expressions for Tangence
    template <class Term>
    struct Tan : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Tan(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::tan(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("tan(") + term<0>(*this).str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for arcus cosine
    template <class Term>
    struct Acos : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Acos(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::acos(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("acos(") + term<0>(*this).str() + ")"; }
    };

    /// Expressions for arcus sine
    template <class Term>
    struct Asin : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Asin(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::asin(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("asin(") + term<0>(*this).str() + ")"; }
    };

    /// Expressions for arcus tangence
    template <class Term>
    struct Atan : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Atan(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::atan(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("atan(") + term<0>(*this).str() + ")"; }
    };

    /// Expressions for arcus tangence2, i.e. atan(x/y)
    template <class Term1, class Term2>
    struct Atan2 : public LazyOperatorTerms<Term1, Term2>
    {
      using Super = LazyOperatorTerms<Term1, Term2>;
      using value_type = Value_t<Term1>;
      int deg;
      
      Atan2(Term1&& term1_, Term2&& term2_, int degree_ = 1) 
        : Super(std::forward<Term1>(term1_), std::forward<Term2>(term2_)), deg(degree_) 
      { }
      
      int getDegree() const
      {
        return deg * (degree<0>(*this) + degree<1>(*this));
      }

      value_type operator()(const int& iq) const { return std::atan2(term<0>(*this)(iq), term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("atan2(") + term<0>(*this).str() + ", " + term<1>(*this) + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for cosine hyperbolicus
    template <class Term>
    struct Cosh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Cosh(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::cosh(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("cosh(") + term<0>(*this).str() + ")"; }
    };

    /// Expressions for sine hyperbolicus
    template <class Term>
    struct Sinh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Sinh(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::sinh(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("sinh(") + term<0>(*this).str() + ")"; }
    };

    /// Expressions for tangence hyperbolicus
    template <class Term>
    struct Tanh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Tanh(Term&& term_, int degree_ = 1) : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * (degree<0>(*this));
      }

      value_type operator()(const int& iq) const { return std::tanh(term<0>(*this)(iq)); }
      
      std::string str() const { return std::string("tanh(") + term<0>(*this).str() + ")"; }
    };
    
    // ___________________________________________________________________________
    

    /// Expressions for arcus cosine hyperbolicus
    template <class Term>
    struct Acosh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Acosh(Term&& term_, int degree_ = 1) 
        : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * degree<0>(*this);
      }

      value_type operator()(const int& iq) const 
      { 
        value_type tmp = term<0>(*this)(iq);
      	return std::log(tmp + std::sqrt(sqr(tmp) - 1.0)); 
      }
      
      std::string str() const { return std::string("acosh(") + term<0>(*this).str() + ")"; }
    };

    /// Expressions for arcus sine hyperbolicus
    template <class Term>
    struct Asinh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Asinh(Term&& term_, int degree_ = 1) 
        : Super(std::forward<Term>(term_)), deg(degree_) {}
      
      int getDegree() const
      {
        return deg * degree<0>(*this);
      }

      value_type operator()(const int& iq) const 
      { 
        value_type tmp = term<0>(*this)(iq);
        return std::log(tmp + std::sqrt(sqr(tmp) + 1.0)); 
      }
      
      std::string str() const { return std::string("asinh(") + term<0>(*this).str() + ")"; }
    };

    /// Expressions for arcus tangence hyperbolicus
    template <class Term>
    struct Atanh : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      int deg;
      
      Atanh(Term&& term_, int deg_ = 1) : Super(std::forward<Term>(term_)), deg(deg_) {}
      
      int getDegree() const
      {
        return deg * degree<0>(*this);
      }

      value_type operator()(const int& iq) const 
      { 
      	value_type tmp = term<0>(*this)(iq);
      	return 0.5 * std::log((1.0 + tmp) / (1.0 - tmp)); 
      }
      
      std::string str() const { return std::string("atanh(") + term<0>(*this).str() + ")"; }
    };
    
    
    // ___________________________________________________________________________
    
    
    /// Expressions for the maximum of two expressions max(E1, E2)
    template <class Term1, class Term2>
    struct Max : public LazyOperatorTerms<Term1, Term2>
    {
      using Super = LazyOperatorTerms<Term1, Term2>;
      using value_type = typename std::common_type<Value_t<Term1>, Value_t<Term2> >::type;
      
      Max(Term1&& term1_, Term2&& term2_)
        : Super(std::forward<Term1>(term1_), std::forward<Term2>(term2_)) 
      { }
      
      int getDegree() const
      {
        return std::max(degree<0>(*this), degree<1>(*this));
      }

      value_type operator()(const int& iq) const
      {
        return std::max(term<0>(*this)(iq), term<1>(*this)(iq)); 
      }
      
      std::string str() const { return std::string("max(") + term<0>(*this).str() + ", " + term<1>(*this).str() + ")"; }
    };
    
    
    /// Expressions for the minimum of two expressions min(E1, E2)
    template <class Term1, class Term2>
    struct Min : public LazyOperatorTerms<Term1, Term2>
    {
      using Super = LazyOperatorTerms<Term1, Term2>;
      using value_type = typename std::common_type<Value_t<Term1>, Value_t<Term2> >::type;
      
      Min(Term1&& term1_, Term2&& term2_)
        : Super(std::forward<Term1>(term1_), std::forward<Term2>(term2_)) 
      { }
      
      int getDegree() const
      {
        return std::max(degree<0>(*this), degree<1>(*this));
      }

      value_type operator()(const int& iq) const 
      {
      	return std::min(term<0>(*this)(iq), term<1>(*this)(iq)); 
      }
      
      std::string str() const { return std::string("min(") + term<0>(*this).str() + ", " + term<1>(*this).str() + ")"; }
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
  template <class Term1, class Term2>
  inline typename result_of::Max<Term1, Term2>::type
  max(Term1&& t1, Term2&& t2) 
  { 
    using Expr1 = traits::to_expr<Term1>;
    using Expr2 = traits::to_expr<Term2>;
    return {Expr1::get(t1), Expr2::get(t2)};
  }

  // minimum of two terms
  // _____________________________________________________________________________
  template <class Term1, class Term2>
  inline typename result_of::Min<Term1, Term2>::type
  min(Term1&& t1, Term2&& t2) 
  { 
    using Expr1 = traits::to_expr<Term1>;
    using Expr2 = traits::to_expr<Term2>;
    return {Expr1::get(t1), Expr2::get(t2)};
  }


  //______________________________________________________________________________

  // absolute value of a term
  template <class Term>
  inline typename boost::enable_if< traits::is_expr<Term>,
    expressions::Abs<Term> >::type
  abs_(Term&& t) { return {std::forward<Term>(t)}; } // TODO: Funktionsnamen ohne Unterstrich

  // signum of a term
  template <class Term>
  inline typename boost::enable_if< traits::is_expr<Term>,
    expressions::Signum<Term> >::type
  signum(Term&& t) { return {std::forward<Term>(t)}; }

  // 
  template <class Term>
  inline typename boost::enable_if< traits::is_expr<Term>,
    expressions::Ceil<Term> >::type
  ceil(Term&& t) { return {std::forward<Term>(t)}; }

  // 
  template <class Term>
  inline typename boost::enable_if< traits::is_expr<Term>,
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
  inline typename boost::enable_if<
    traits::is_expr<Term>,
    expressions::Atanh<Term> >::type
  atanh(Term&& t) { return {std::forward<Term>(t)}; }

} // end namespace AMDiS
